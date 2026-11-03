"""Annotation report: figures + tables + one self-contained HTML page.

Consumes the outputs of ``pcip.predict.predict_og_adata`` (per-cell predictions and per-type scores) together
with the OG-layer AnnData that was scored and the dictionary used. Produces, in ``out_dir``:

  {id}.composition.tsv      predicted cell-type counts / fractions, median margin, tie and invalid counts
  {id}.markers.tsv          per type: its program OGs (weight, query member genes, detection in-type vs rest,
                            mean expression in-type vs rest, log2 fold change): the marker list
  {id}.group_scores.tsv     mean per-type program score within each predicted group (the score matrix)
  {id}.umap.tsv             2-D embedding coordinates per cell
  figures/*.png, *.pdf      UMAPs by predicted type, tissue system and confidence, composition, score heatmap,
                            marker dot plot (query gene ids), margin per type, state scores
  {id}.report.html          everything above in one page (figures embedded; no external files needed)

The embedding is computed on the OG matrix (normalise + log1p unless the input is already scaled, top-variance
features, PCA, kNN, UMAP; fixed seed), so it is comparable across species. Nothing here changes predictions.
"""
from __future__ import annotations

import base64
import html
import json
from pathlib import Path
from typing import Any

import numpy as np

N_MARKERS_PLOT = 5
N_HVG = 2000


def _dense(X) -> np.ndarray:
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def _is_scaled(X) -> bool:
    m = X.min() if hasattr(X, "min") else np.min(X)
    return bool(np.asarray(m).item() < 0)


def top_variance_features(X, n: int) -> np.ndarray:
    """Indices of the n highest-variance columns; sparse input is never densified as a whole."""
    if hasattr(X, "toarray"):
        mean = np.asarray(X.mean(axis=0)).ravel()
        var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean ** 2
    else:
        var = np.asarray(X).var(axis=0)
    return np.argsort(-var, kind="stable")[:min(n, X.shape[1])]


def embed(og_adata: Any, seed: int, n_pcs: int = 30) -> np.ndarray:
    """2-D UMAP of the OG matrix (seeded)."""
    import scanpy as sc
    a = og_adata.copy()
    if not _is_scaled(a.X):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    a = a[:, top_variance_features(a.X, N_HVG)].copy()
    sc.pp.pca(a, n_comps=max(2, min(n_pcs, a.n_vars - 1, a.n_obs - 1)), random_state=seed)
    sc.pp.neighbors(a, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    return np.asarray(a.obsm["X_umap"])


def composition(preds) -> "Any":
    import pandas as pd
    p = preds["predicted_cell_type"].astype(str)
    rows = []
    for t, g in preds.groupby(p, sort=True):
        rows.append({"cell_type": t if t else "(invalid)", "n_cells": len(g), "fraction": len(g) / len(preds),
                     "median_margin": float(np.nanmedian(g["margin"])) if len(g) else np.nan,
                     "n_tie_at_top": int(g["tie_at_top"].sum()), "n_invalid": int(g["invalid_score"].sum())})
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False).reset_index(drop=True)


def marker_table(og_adata: Any, preds, dictionary: dict) -> "Any":
    """Per type, every program OG with weight, query member genes and in-type vs rest expression stats."""
    import pandas as pd
    X = og_adata.X
    lib = np.asarray(X.sum(axis=1)).ravel()
    scaled = _is_scaled(X)
    var_index = {str(v): i for i, v in enumerate(og_adata.var_names)}
    members = (og_adata.var["member_gene_ids"].astype(str).to_numpy() if "member_gene_ids" in og_adata.var
               else np.array([""] * og_adata.n_vars))
    pred = preds["predicted_cell_type"].astype(str).to_numpy()
    rows = []
    for p in dictionary.get("programs", []):
        t = str(p.get("cell_type_coarse") or p.get("cell_type") or p.get("program_id"))
        inn = pred == t
        out = ~inn
        for rank, f in enumerate(p.get("positive_features", []), start=1):
            fid = str(f["feature_id"])
            j = var_index.get(fid)
            row = {"cell_type": t, "program_id": p.get("program_id"), "rank": rank, "orthogroup_id": fid,
                   "weight": float(f.get("weight", 1.0)), "present_in_query": j is not None,
                   "query_member_genes": members[j] if j is not None else ""}
            if j is not None:
                col = _dense(X[:, j]).ravel()
                both = inn.any() and out.any()
                if scaled:   # centred values: "detected" and fold change are undefined; report means only
                    row.update({"detect_frac_in_type": np.nan, "detect_frac_rest": np.nan,
                                "mean_expr_in_type": float(col[inn].mean()) if inn.any() else np.nan,
                                "mean_expr_rest": float(col[out].mean()) if out.any() else np.nan,
                                "log2fc_in_vs_rest": np.nan})
                else:
                    norm = col / np.maximum(lib, 1) * 1e4          # per-cell normalised expression (CP10k)
                    expr = np.log1p(norm)
                    row.update({
                        "detect_frac_in_type": float((col[inn] > 0).mean()) if inn.any() else np.nan,
                        "detect_frac_rest": float((col[out] > 0).mean()) if out.any() else np.nan,
                        "mean_expr_in_type": float(expr[inn].mean()) if inn.any() else np.nan,
                        "mean_expr_rest": float(expr[out].mean()) if out.any() else np.nan,
                        # log2 ratio of group mean normalised expression (not of back-transformed mean logs)
                        "log2fc_in_vs_rest": float(np.log2((norm[inn].mean() + 1e-3) / (norm[out].mean() + 1e-3)))
                                             if both else np.nan})
            rows.append(row)
    return pd.DataFrame(rows)


def group_scores(preds, scores) -> "Any":
    g = preds["predicted_cell_type"].astype(str).to_numpy()
    df = scores.copy()
    df.columns = [c.replace("score_", "") for c in df.columns]
    df["predicted"] = g
    return df.groupby("predicted").mean(numeric_only=True)


# Fixed colours so a cell type looks the same in every run and organ; related types share a hue family.
TYPE_COLORS = {
    "epidermis": "#6aa84f", "cortex": "#e3a23b", "endodermis": "#a8612a", "pericycle": "#8e6bb8",
    "procambium": "#69b3d8", "phloem": "#d9577b", "xylem": "#2f5d9e", "root_cap": "#c7b356",
    "mesophyll": "#3d8b40", "vascular": "#2f5d9e",
    "pollen": "#e9c43f", "stamen": "#e0853a", "carpel": "#b8395d", "perianth": "#d88fbf", "non_floral": "#9a9a9a",
    "embryo": "#4c9a6a", "endosperm": "#e6c06b", "seed_coat": "#8b5a3c",
    "dermal": "#6aa84f", "ground": "#e3a23b",
}
CONF_COLORS = ["#d9d9d9", "#9ecae1", "#4292c6", "#08306b"]
FONT_FAMILY = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]


def _palette(types):
    import matplotlib
    cmap = matplotlib.colormaps.get_cmap("tab20")
    extra = [t for t in sorted(types) if t not in TYPE_COLORS]
    pal = {t: TYPE_COLORS[t] for t in types if t in TYPE_COLORS}
    pal.update({t: cmap(i % 20) for i, t in enumerate(extra)})
    return pal


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": FONT_FAMILY, "font.size": 9,
                         "axes.titlesize": 10, "axes.titleweight": "bold", "axes.labelsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.fontsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    return plt


def _save(plt, fig, fig_dir: Path, key: str, out: dict) -> None:
    out[key] = fig_dir / f"{key}.png"
    fig.savefig(out[key], dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / f"{key}.pdf", bbox_inches="tight")
    plt.close(fig)


def _umap_axes(ax) -> None:
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _categorical_umap(plt, xy, labels, pal):
    labels = np.asarray(labels).astype(str)
    cats = [c for c in sorted(set(labels) - {""}, key=lambda c: -(labels == c).sum())]
    fig, ax = plt.subplots(figsize=(6.2, 5))
    order = np.random.default_rng(0).permutation(len(labels))   # mixed draw order so no type hides another
    ax.scatter(xy[order, 0], xy[order, 1], s=3, c=[pal.get(labels[i], "#cccccc") for i in order],
               linewidths=0, rasterized=True)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=6, color=pal[c],
                          label=f"{c.replace('_', ' ')}  ({(labels == c).sum():,})") for c in cats]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, handletextpad=0.3)
    _umap_axes(ax)
    return fig


def _marker_label(row) -> str:
    genes = [g for g in str(row["query_member_genes"]).split(",") if g and g != "nan"]
    if not genes:
        return str(row["orthogroup_id"]).split(":")[-1]
    return genes[0] + (f" +{len(genes) - 1}" if len(genes) > 1 else "")


def make_figures(fig_dir: Path, xy, preds, comp, gscores, markers, og_adata,
                 reliability: dict | None = None) -> dict[str, Path]:
    plt = _style()
    fig_dir.mkdir(parents=True, exist_ok=True)
    pred = preds["predicted_cell_type"].astype(str).to_numpy()
    types = [t for t in comp["cell_type"] if t != "(invalid)"]          # most frequent first
    pal = _palette(types)
    out: dict[str, Path] = {}

    fig = _categorical_umap(plt, xy, pred, pal)
    _save(plt, fig, fig_dir, "umap_predicted", out)

    if "tissue_system" in preds:
        ts = preds["tissue_system"].astype(str).to_numpy()
        fig = _categorical_umap(plt, xy, ts, _palette(sorted(set(ts) - {""})))
        _save(plt, fig, fig_dir, "umap_tissue_system", out)

    acc = (reliability or {}).get("margin_quartile_accuracy")
    if "confidence_quartile" in preds:
        q = preds["confidence_quartile"].to_numpy(dtype=int)
        fig, ax = plt.subplots(figsize=(6.2, 5))
        for k in (1, 2, 3, 4):
            m = q == k
            lab = f"Q{k}" + (f"  (expected accuracy {acc[k - 1]:.0%})" if acc else "") + f"  ({m.sum():,})"
            ax.scatter(xy[m, 0], xy[m, 1], s=3, color=CONF_COLORS[k - 1], linewidths=0, rasterized=True, label=lab)
        ax.legend(title="confidence quartile", markerscale=3, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  frameon=False, title_fontsize=8)
    else:
        fig, ax = plt.subplots(figsize=(6.2, 5))
        sca = ax.scatter(xy[:, 0], xy[:, 1], s=3, c=preds["margin"].to_numpy(), cmap="viridis", linewidths=0,
                         rasterized=True)
        fig.colorbar(sca, ax=ax, label="margin (top - runner-up score)", shrink=0.7)
    _umap_axes(ax)
    _save(plt, fig, fig_dir, "umap_margin", out)

    # composition, split by confidence quartile when available
    fig, ax = plt.subplots(figsize=(6.2, 0.32 * len(types) + 1.2))
    yy = np.arange(len(types))[::-1]
    n = np.array([(pred == t).sum() for t in types])
    if "confidence_quartile" in preds:
        q = preds["confidence_quartile"].to_numpy(dtype=int)
        left = np.zeros(len(types))
        for k in (1, 2, 3, 4):
            w = np.array([((pred == t) & (q == k)).sum() for t in types]) / len(pred)
            ax.barh(yy, w, left=left, color=CONF_COLORS[k - 1], height=0.7, label=f"Q{k}")
            left += w
        ax.legend(title="confidence quartile (Q1 = least confident)", ncol=4, loc="upper center", frameon=False,
                  title_fontsize=8, bbox_to_anchor=(0.5, -0.10 - 0.5 / len(types)))
    else:
        ax.barh(yy, n / len(pred), color=[pal[t] for t in types], height=0.7)
    for y, k in zip(yy, n):
        ax.text(k / len(pred), y, f"  {k:,} ({k / len(pred):.0%})", va="center", fontsize=8)
    ax.set_yticks(yy); ax.set_yticklabels([t.replace("_", " ") for t in types])
    for lab, t in zip(ax.get_yticklabels(), types):
        lab.set_color(pal[t]); lab.set_fontweight("bold")
    ax.set_xlabel("fraction of cells"); ax.set_xlim(0, max(n / len(pred)) * 1.3)
    _save(plt, fig, fig_dir, "composition", out)

    # mean program score per predicted group (scores are centred on the dataset background, so 0 = background)
    G = gscores.loc[[t for t in types if t in gscores.index], [t for t in types if t in gscores.columns]
                    + [c for c in gscores.columns if c not in types]]
    vals = G.to_numpy(dtype=float)
    lim = float(np.nanmax(np.abs(vals))) or 1.0
    fig, ax = plt.subplots(figsize=(0.5 * G.shape[1] + 2.4, 0.42 * G.shape[0] + 1.6))
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    if G.size <= 400:
        for i in range(G.shape[0]):
            for j in range(G.shape[1]):
                ax.text(j, i, f"{vals[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(vals[i, j]) > 0.6 * lim else "#333333")
    ax.set_yticks(range(G.shape[0])); ax.set_yticklabels([t.replace("_", " ") for t in G.index])
    ax.set_xticks(range(G.shape[1])); ax.set_xticklabels([t.replace("_", " ") for t in G.columns], rotation=45,
                                                          ha="right")
    ax.set_ylabel("predicted cell type"); ax.set_xlabel("program score")
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    fig.colorbar(im, ax=ax, label="mean score (0 = dataset background)", shrink=0.8)
    _save(plt, fig, fig_dir, "score_heatmap", out)

    # top program genes per type, labelled with the query's own gene ids
    mk = markers[markers["present_in_query"]].copy()
    if len(mk) and "detect_frac_in_type" in mk and mk["detect_frac_in_type"].notna().any():
        mk = mk[(mk["detect_frac_in_type"].fillna(0) + mk["detect_frac_rest"].fillna(0)) > 0]
    mk = (mk[mk["cell_type"].isin(types)].assign(_o=lambda d: d["cell_type"].map({t: i for i, t in enumerate(types)}))
          .sort_values(["_o", "rank"]).drop_duplicates("orthogroup_id").groupby("cell_type", sort=False)
          .head(N_MARKERS_PLOT))
    if len(mk) and "mean_expr_in_type" in mk:
        feats = list(mk["orthogroup_id"])
        var_index = {str(v): i for i, v in enumerate(og_adata.var_names)}
        X = og_adata.X
        lib = np.asarray(X.sum(axis=1)).ravel()
        scaled = _is_scaled(X)
        frac = np.zeros((len(types), len(feats))); mean = np.zeros_like(frac)
        for j, f in enumerate(feats):
            col = _dense(X[:, var_index[f]]).ravel()
            expr = col if scaled else np.log1p(col / np.maximum(lib, 1) * 1e4)
            for i, t in enumerate(types):
                m = pred == t
                frac[i, j] = (0.5 if scaled else (col[m] > 0).mean()) if m.any() else 0
                mean[i, j] = expr[m].mean() if m.any() else 0
        rng = mean.max(axis=0) - mean.min(axis=0)
        mscaled = (mean - mean.min(axis=0)) / np.where(rng > 0, rng, 1)   # per gene, 0-1 across types
        fig, ax = plt.subplots(figsize=(max(5.0, 0.24 * len(feats) + 2.6), 0.34 * len(types) + 2.4))
        yy, xx = np.meshgrid(range(len(types)), range(len(feats)), indexing="ij")
        sca = ax.scatter(xx.ravel(), yy.ravel(), s=frac.ravel() * 90 + 0.5, c=mscaled.ravel(), cmap="Reds",
                         vmin=0, vmax=1, edgecolors="#555555", linewidths=0.2)
        ax.set_xticks(range(len(feats))); ax.set_xticklabels([_marker_label(r) for _, r in mk.iterrows()],
                                                              rotation=90, fontsize=6.5)
        ax.set_yticks(range(len(types))); ax.set_yticklabels([t.replace("_", " ") for t in types])
        ax.set_ylim(len(types) - 0.5, -1.9); ax.set_xlim(-0.6, len(feats) - 0.4)
        start = 0
        for t, g in mk.groupby("cell_type", sort=False):
            k = len(g)
            ax.plot([start - 0.35, start + k - 0.65], [-1.0, -1.0], color=pal[t], lw=4, solid_capstyle="butt",
                    clip_on=False)
            ax.text(start + (k - 1) / 2, -1.25, t.replace("_", " "), ha="center", va="bottom", fontsize=7.5,
                    color=pal[t], fontweight="bold")
            start += k
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        cb = fig.colorbar(sca, ax=ax, shrink=0.4, pad=0.02, anchor=(0.0, 1.0))
        cb.set_label("mean expression\n(scaled per gene)")
        if not scaled:
            for v in (0.25, 0.5, 0.75, 1.0):
                ax.scatter([], [], s=v * 90 + 0.5, color="#bbbbbb", edgecolors="#555555", linewidths=0.2,
                           label=f"{v:.0%}")
            ax.legend(title="cells expressing", loc="upper left", bbox_to_anchor=(1.01, 0.42), frameon=False,
                      title_fontsize=7, fontsize=7, labelspacing=0.9, borderaxespad=0)
        ax.set_xlabel("query genes of the top program orthogroups (coloured bar = program's cell type)"
                      + ("\nscaled input: dot size not meaningful" if scaled else ""))
        _save(plt, fig, fig_dir, "marker_dotplot", out)

    fig, ax = plt.subplots(figsize=(6.2, 0.32 * len(types) + 1.2))
    data = [preds.loc[pred == t, "margin"].dropna().to_numpy() for t in types][::-1]
    bp = ax.boxplot(data, vert=False, showfliers=False, patch_artist=True, widths=0.6,
                    medianprops={"color": "#222222"})
    for patch, t in zip(bp["boxes"], types[::-1]):
        patch.set_facecolor(pal[t]); patch.set_alpha(0.85)
    ax.set_yticks(range(1, len(types) + 1)); ax.set_yticklabels([t.replace("_", " ") for t in types[::-1]])
    ax.axvline(0, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("margin (top - runner-up score)")
    _save(plt, fig, fig_dir, "margin_by_type", out)
    return out


SAMPLE_KEYS = ("sample", "orig.ident", "batch", "library", "replicate")
CLUSTER_KEYS = ("seurat_clusters", "leiden", "louvain", "clusters", "cluster")
MIN_CELLS_PER_SAMPLE = 50
# Mean total-variation distance between the samples' call distributions within one cluster. On replicate samples of
# labelled root atlases (soybean, Arabidopsis, rice) it was 0.04-0.14; on a dataset whose calls disagreed between
# biological replicates it was 0.27-0.58.
TV_CONSISTENT, TV_UNRELIABLE = 0.20, 0.30
MIN_ASSESSED_FRACTION = 0.20      # below this share of cells in comparable clusters the check is not reported
MAX_CLUSTER_FRACTION = 0.50       # a cluster column whose largest cluster holds more cells than this is ignored
# Share of (cluster, other sample) pairs whose expression profile is closest to the same cluster in the other sample.
# Integrated clusters of labelled root atlases: 0.99-1.00; clusters that did not align across replicates: 0.12.
MIN_CLUSTER_CORRESPONDENCE = 0.80


def pick_obs_key(obs, given: str | None, candidates: tuple[str, ...], max_levels: int = 500) -> str | None:
    """`given` if set (must exist), else the first candidate column with 2..max_levels values."""
    if given:
        if given not in obs.columns:
            raise SystemExit(f"report: obs column {given!r} not found (columns: {', '.join(map(str, obs.columns))})")
        return given
    for c in candidates:
        if c in obs.columns and 2 <= obs[c].nunique() <= max_levels:
            return c
    return None


def _tv(a, b) -> float:
    idx = a.index.union(b.index)
    return 0.5 * float((a.reindex(idx, fill_value=0) - b.reindex(idx, fill_value=0)).abs().sum())


def usable_cluster_key(obs, key: str | None) -> str | None:
    """Drop a cluster column that does not partition the cells (e.g. stale metadata with one dominant cluster)."""
    if key is None:
        return None
    return key if obs[key].astype(str).value_counts(normalize=True).iloc[0] <= MAX_CLUSTER_FRACTION else None


def cluster_correspondence(X, samples, clusters, n_features: int = N_HVG) -> float | None:
    """Do the clusters mean the same cells in every sample? Pseudobulk each (sample, cluster) with
    >= MIN_CELLS_PER_SAMPLE cells; for every cluster and every other sample, check whether the most correlated
    cluster of that sample is the same cluster. Returns the cell-weighted share of hits (None if nothing to compare)."""
    import itertools
    import pandas as pd
    import scipy.sparse as sp
    key = pd.MultiIndex.from_arrays([np.asarray(samples).astype(str), np.asarray(clusters).astype(str)])
    codes, groups = pd.factorize(key)
    n = np.bincount(codes)
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(groups), len(codes)))
    pb = _dense(G @ X)
    if _is_scaled(X):
        prof = pb / n[:, None]
    else:
        prof = np.log1p(pb / np.maximum(pb.sum(1, keepdims=True), 1) * 1e6)
    keep = n >= MIN_CELLS_PER_SAMPLE
    prof, n = prof[keep], n[keep]
    smp = np.asarray(groups.get_level_values(0))[keep]
    clu = np.asarray(groups.get_level_values(1))[keep]
    prof = prof[:, np.argsort(-prof.var(axis=0), kind="stable")[:n_features]]
    C = np.corrcoef(prof)
    hit = total = 0.0
    for s, t in itertools.permutations(np.unique(smp), 2):
        it = np.flatnonzero(smp == t)
        if len(it) < 2:
            continue
        for i in np.flatnonzero(smp == s):
            if clu[i] not in clu[it]:
                continue
            hit += n[i] * (clu[it[np.argmax(C[i, it])]] == clu[i])
            total += n[i]
    return float(hit / total) if total else None


def sample_composition_tv(calls, samples) -> float | None:
    """Median pairwise total-variation distance between the samples' overall call distributions."""
    import itertools
    import pandas as pd
    d = pd.DataFrame({"call": np.asarray(calls).astype(str), "sample": np.asarray(samples).astype(str)})
    ok = [s for s, n in d["sample"].value_counts().items() if n >= MIN_CELLS_PER_SAMPLE]
    dist = {s: d.loc[d["sample"] == s, "call"].value_counts(normalize=True) for s in ok}
    tvs = [_tv(dist[a], dist[b]) for a, b in itertools.combinations(ok, 2)]
    return float(np.median(tvs)) if tvs else None


def sample_consistency(calls, samples, clusters) -> tuple[Any, dict]:
    """Do cells of one cluster get the same calls in every sample?

    For each cluster with >= MIN_CELLS_PER_SAMPLE cells in at least two samples: the call distribution per sample,
    their mean pairwise total-variation distance (0 = identical, 1 = disjoint) and whether the most frequent call is
    the same in every sample. Clusters, not author labels, so it runs on any dataset."""
    import itertools
    import pandas as pd
    d = pd.DataFrame({"call": np.asarray(calls).astype(str), "sample": np.asarray(samples).astype(str),
                      "cluster": np.asarray(clusters).astype(str)})
    rows = []
    for cl, g in d.groupby("cluster", sort=False):
        ok = [s for s, n in g["sample"].value_counts().items() if n >= MIN_CELLS_PER_SAMPLE]
        if len(ok) < 2:
            continue
        dist = {s: g.loc[g["sample"] == s, "call"].value_counts(normalize=True) for s in ok}
        tops = sorted({v.index[0] for v in dist.values()})
        rows.append({"cluster": cl, "n_cells": len(g), "n_samples": len(ok),
                     "mean_tv": float(np.mean([_tv(dist[a], dist[b]) for a, b in itertools.combinations(ok, 2)])),
                     "same_top_call": len(tops) == 1, "top_calls": "/".join(tops)})
    table = pd.DataFrame(rows, columns=["cluster", "n_cells", "n_samples", "mean_tv", "same_top_call", "top_calls"])
    if table.empty:
        return table, {"status": "not_assessable",
                       "reason": f"no cluster has >= {MIN_CELLS_PER_SAMPLE} cells in two or more samples"}
    w = table["n_cells"].to_numpy(float)
    assessed = float(w.sum() / len(d))
    if assessed < MIN_ASSESSED_FRACTION:
        return table, {"status": "not_assessable", "fraction_cells_assessed": round(assessed, 4),
                       "reason": f"only {assessed:.0%} of the cells are in clusters shared by two or more samples "
                                 "(clusters are sample-specific; integrate the samples and pass --cluster-key)"}
    tv = float(np.average(table["mean_tv"], weights=w))
    status = "consistent" if tv <= TV_CONSISTENT else ("check" if tv <= TV_UNRELIABLE else "not_reproducible")
    return table.sort_values("n_cells", ascending=False).reset_index(drop=True), {
        "status": status, "weighted_mean_tv": round(tv, 4), "fraction_cells_assessed": round(assessed, 4),
        "cells_in_clusters_with_same_top_call": round(float(w[table["same_top_call"].to_numpy()].sum() / w.sum()), 4),
        "n_clusters_tested": int(len(table)), "thresholds": {"consistent": TV_CONSISTENT, "unreliable": TV_UNRELIABLE}}


def margin_quartile(margin) -> "Any":
    """Per-cell confidence: quartile 1-4 of the cell's margin within this query (4 = largest margins). On held-out
    species, accuracy rose monotonically across these quartiles, whereas per-type trust levels did not transfer."""
    m = np.asarray(margin, dtype=float)
    finite = np.isfinite(m)
    if not finite.any():
        return np.ones(m.shape, dtype=int)
    qs = np.quantile(m[finite], [0.25, 0.5, 0.75])
    return np.where(finite, np.digitize(np.nan_to_num(m, nan=-np.inf), qs) + 1, 1).astype(int)


def confidence_table(preds, reliability: dict | None) -> "Any":
    """Cells per margin quartile + the accuracy measured for that quartile on held-out species
    (``reliability.margin_quartile_accuracy``; NaN when the dictionary carries no measurement)."""
    import pandas as pd
    q = margin_quartile(preds["margin"])
    acc = (reliability or {}).get("margin_quartile_accuracy") or [float("nan")] * 4
    return pd.DataFrame({"confidence_quartile": [1, 2, 3, 4],
                         "n_cells": [int((q == k).sum()) for k in (1, 2, 3, 4)],
                         "expected_accuracy": [acc[k - 1] for k in (1, 2, 3, 4)]})


def reliability_table(comp, reliability: dict) -> "Any":
    """Per predicted type: cell count + cross-species (leave-one-species-out) F1 from the dictionary's training-time
    evaluation. Reference only: a type's F1 on other species did not predict its accuracy on a new one."""
    import pandas as pd
    f1 = (reliability or {}).get("coarse_loso_f1", {})
    rows = [{"cell_type": r["cell_type"], "n_cells": int(r["n_cells"]), "fraction": r["fraction"],
             "cross_species_f1": f1.get(r["cell_type"])}
            for _, r in comp[comp["cell_type"] != "(invalid)"].iterrows()]
    return pd.DataFrame(rows, columns=["cell_type", "n_cells", "fraction", "cross_species_f1"])


def make_extra_figures(fig_dir: Path, xy, full_preds) -> dict[str, Path]:
    """One UMAP per developmental-state score (fine types are not shown)."""
    plt = _style()
    out: dict[str, Path] = {}
    for c in [c for c in full_preds.columns if c.startswith("state_") and c.endswith("_score")]:
        name = c[len("state_"):-len("_score")]
        label = name.split("_qc")[0].replace("_", " ")
        v = full_preds[c].to_numpy(dtype=float)
        order = np.argsort(np.nan_to_num(v, nan=-np.inf))
        fig, ax = plt.subplots(figsize=(6.2, 5))
        sca = ax.scatter(xy[order, 0], xy[order, 1], s=3, c=v[order], cmap="viridis", linewidths=0, rasterized=True)
        fig.colorbar(sca, ax=ax, label=f"{label} state score", shrink=0.7)
        _umap_axes(ax)
        _save(plt, fig, fig_dir, f"umap_state_{name}", out)
    return out


def _img(path: Path, caption: str = "") -> str:
    data = base64.b64encode(path.read_bytes()).decode()
    cap = f"<figcaption>{caption}</figcaption>" if caption else ""
    return (f'<figure><a href="figures/{path.name}"><img src="data:image/png;base64,{data}" alt="{path.stem}"></a>'
            f"{cap}</figure>")


def _table(df, n=None) -> str:
    d = df if n is None else df.head(n)
    return d.to_html(index=False, float_format=lambda x: f"{x:.3g}", border=0, classes="t")


_CSS = """
body{font-family:Arial,Helvetica,sans-serif;max-width:1180px;margin:auto;padding:24px;color:#222;line-height:1.45}
h1{font-size:24px;margin-bottom:4px}h2{font-size:18px;margin-top:36px;border-bottom:2px solid #e5e5e5;padding-bottom:4px}
.sub{color:#666;font-size:13px}.cards{display:flex;flex-wrap:wrap;gap:12px;margin:18px 0}
.card{flex:1 1 150px;background:#f6f8f6;border-radius:8px;padding:10px 14px}
.card b{display:block;font-size:22px;color:#2f5d33}.card span{font-size:12px;color:#555}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:18px}
figure{margin:0}figure img{width:100%;border:1px solid #eee;border-radius:4px}
figcaption{font-size:12px;color:#555;margin-top:4px}
table.t{border-collapse:collapse;font-size:12px;margin:8px 0}table.t td,table.t th{padding:3px 9px;border-bottom:1px solid #e5e5e5;text-align:left}
table.t th{background:#f3f3f3}.note{background:#fff8e6;border-left:4px solid #e3a23b;padding:8px 12px;font-size:13px}
.alert{background:#fdecea;border-left:4px solid #c62828;padding:10px 14px;font-size:14px;margin:12px 0}
.meta{font-size:11px;color:#555;white-space:pre-wrap;background:#fafafa;padding:8px}code{background:#f3f3f3;padding:0 3px}
"""


def write_html(path: Path, title: str, meta: dict, comp, markers, figs: dict[str, Path], rel=None, reliability=None,
               conf=None, consistency: dict | None = None, cons_table=None) -> None:
    top = (markers[markers["present_in_query"]].sort_values(["cell_type", "rank"]).groupby("cell_type").head(10)
           [["cell_type", "rank", "orthogroup_id", "weight", "query_member_genes", "detect_frac_in_type",
             "detect_frac_rest", "log2fc_in_vs_rest"]] if "detect_frac_in_type" in markers else markers.head(0))
    n_cells = int(comp["n_cells"].sum())
    acc = (reliability or {}).get("margin_quartile_accuracy")
    cards = [(f"{n_cells:,}", "cells annotated"), (str(int((comp["cell_type"] != "(invalid)").sum())), "cell types found"),
             (f"{comp['cell_type'].iloc[0].replace('_', ' ')}", f"largest type ({comp['fraction'].iloc[0]:.0%})")]
    if "species" in meta:
        cards.append((html.escape(str(meta["species"])), "species"))
    if acc:
        cards.append((f"{acc[-1]:.0%} / {acc[0]:.0%}", "expected accuracy, most / least confident quarter"))
    card_html = "".join(f"<div class='card'><b>{v}</b><span>{k}</span></div>" for v, k in cards)

    cons = consistency or {"status": "not_run"}
    banner = ""
    if cons["status"] == "not_reproducible":
        banner = ("<div class='alert'><b>The calls are not reproducible across samples.</b> Cells of the same cluster "
                  f"received different cell types in different samples (mean difference {cons['weighted_mean_tv']:.2f}; "
                  f"replicates of labelled root datasets: 0.04–0.14). Do not use these calls without independent "
                  "confirmation. See <i>Reproducibility across samples</i>.</div>")
    elif cons["status"] == "check":
        banner = ("<div class='note'>Calls differ noticeably between samples within the same cluster "
                  f"(mean difference {cons['weighted_mean_tv']:.2f}). Check the per-cluster table below.</div>")
    if cons["status"] in ("consistent", "check", "not_reproducible"):
        cons_html = ["<h2>Reproducibility across samples</h2>",
                     f"<p>Samples: <code>{html.escape(str(cons['sample_key']))}</code>; clusters: "
                     f"<code>{html.escape(str(cons['cluster_key']))}</code>. For every cluster with at least "
                     f"{MIN_CELLS_PER_SAMPLE} cells in two or more samples, the distribution of calls is compared between "
                     "samples (total-variation distance: 0 = identical, 1 = no overlap). Weighted mean "
                     f"<b>{cons['weighted_mean_tv']:.2f}</b> over {cons['fraction_cells_assessed']:.0%} of the cells; "
                     f"{cons['cells_in_clusters_with_same_top_call']:.0%} of the assessed "
                     "cells are in clusters whose most frequent call is the same in every sample. ≤ "
                     f"{TV_CONSISTENT:.2f} = consistent, > {TV_UNRELIABLE:.2f} = not reproducible. Differences between "
                     "biologically different samples (e.g. time points) can be real; differences between replicates "
                     "are not. The clusters must come from an integration of the samples.</p>", _table(cons_table, 20)]
    else:
        cons_html = ["<h2>Reproducibility across samples</h2>",
                     f"<p>Not assessed: {html.escape(str(cons.get('reason', '')))}.</p>"]
    if cons.get("median_sample_composition_tv") is not None:
        cons_html.append(f"<p>Overall composition: median difference between samples "
                         f"<b>{cons['median_sample_composition_tv']:.2f}</b> (same scale). Samples from different "
                         "conditions or time points can differ for biological reasons.</p>")

    def fig(key, caption):
        return _img(figs[key], caption) if key in figs else ""

    parts = [f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>"
             f"<style>{_CSS}</style></head><body>",
             f"<h1>{html.escape(title)}</h1>",
             "<div class='sub'>One predicted cell type per cell from cross-species gene-program scores. "
             "Every cell receives a call; use the confidence quartile to decide which calls to trust.</div>",
             banner,
             f"<div class='cards'>{card_html}</div>",
             "<h2>Cell types</h2><div class='grid'>",
             fig("umap_predicted", "Predicted cell type of every cell on a UMAP of the orthogroup expression."),
             fig("umap_tissue_system", "The same cells coloured by tissue system."),
             "</div>",
             "<h2>Composition and confidence</h2><div class='grid'>",
             fig("composition", "Fraction of cells per predicted type; shading = confidence quartile."),
             fig("umap_margin", "Prediction confidence: quartile of each cell's margin between its top two program scores "
                               "(Q4 = most confident), with the accuracy measured for that quartile on held-out species."),
             "</div>",
             *(["<p>Cells are ranked by <b>margin</b> (top score − runner-up score) within this dataset and split "
                "into quartiles. "
                + ("On species the dictionary had not seen, calls in quartiles 1–4 matched the author label "
                   + " / ".join(f"{v:.0%}" for v in acc) + " of the time. " if acc else "")
                + "Check calls in quartiles 1–2 with known markers before relying on them.</p>",
                _table(conf)] if conf is not None else []),
             _table(comp),
             "<h2>Program scores</h2><div class='grid'>",
             fig("score_heatmap", "Mean program score of each predicted type; a clear diagonal means the types are "
                                  "well separated."),
             fig("margin_by_type", "Confidence per predicted type: margin (top − runner-up score) of the cells "
                                   "assigned to each type; boxes = quartiles, whiskers = 1.5 × IQR."),
             "</div>",
             *(["<h2>Developmental state</h2><div class='grid'>"]
               + [fig(k, "State score of each cell. It ranks cells and is not a cell type.")
                  for k in figs if k.startswith("umap_state_")] + ["</div>"]
               if any(k.startswith("umap_state_") for k in figs) else []),
             *cons_html,
             "<h2>Marker genes</h2>",
             fig("marker_dotplot", "Your genes in the top-weighted orthogroups of each type's program. Dot size = "
                                   "fraction of cells expressing; colour = mean expression scaled per gene."),
             _table(top), "<p>Full ranked list: <code>markers.tsv</code>.</p>",
             *(["<h2>Reference: cross-species F1 per type</h2>",
                "<div class='note'>F1 of each type on species the dictionary had not seen. It did not reliably "
                "predict how well a type is called in a new species, so it is shown for reference only.</div>",
                _table(rel),
                (f"<p>State scores (cross-species AUROC: "
                 + ", ".join(f"{k.split('_qc')[0]} {v:.2f}" for k, v in reliability["state_auroc"].items()) + ").</p>"
                 if reliability and reliability.get("state_auroc") else "")] if rel is not None else []),
             "<h2>Run details</h2>", f"<div class='meta'>{html.escape(json.dumps(meta, indent=1, default=str))}</div>",
             "</body></html>"]
    path.write_text("\n".join(parts), encoding="utf-8")


def build_report(og_adata: Any, preds, scores, dictionary: dict, out_dir: str | Path, dataset_id: str,
                 seed: int, meta: dict | None = None, full_preds=None, sample_key: str | None = None,
                 cluster_key: str | None = None) -> dict[str, str]:
    """`preds` drive the coarse-level figures / markers; `full_preds` (optional, hierarchical output with
    state_*_score columns) adds the developmental-state section. Fine types are never reported."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    comp = composition(preds)
    markers = marker_table(og_adata, preds, dictionary)
    gsc = group_scores(preds, scores)
    sample_key = pick_obs_key(og_adata.obs, sample_key, SAMPLE_KEYS)
    given_cluster_key = cluster_key
    cluster_key = pick_obs_key(og_adata.obs, cluster_key, CLUSTER_KEYS) if sample_key else None
    if cluster_key and not given_cluster_key:
        cluster_key = usable_cluster_key(og_adata.obs, cluster_key)
    xy = embed(og_adata, seed)
    paths = {"composition": out_dir / f"{dataset_id}.composition.tsv",
             "markers": out_dir / f"{dataset_id}.markers.tsv",
             "group_scores": out_dir / f"{dataset_id}.group_scores.tsv",
             "umap": out_dir / f"{dataset_id}.umap.tsv",
             "html": out_dir / f"{dataset_id}.report.html"}
    comp.to_csv(paths["composition"], sep="\t", index=False)
    markers.to_csv(paths["markers"], sep="\t", index=False)
    gsc.to_csv(paths["group_scores"], sep="\t")
    pd.DataFrame({"cell_id": preds["cell_id"], "umap1": xy[:, 0], "umap2": xy[:, 1]}).to_csv(paths["umap"], sep="\t", index=False)
    reliability = dictionary.get("reliability")
    figs = make_figures(out_dir / "figures", xy, full_preds if full_preds is not None else preds, comp, gsc,
                        markers, og_adata, reliability=reliability)
    if full_preds is not None:
        figs.update(make_extra_figures(out_dir / "figures", xy, full_preds))
    rel = reliability_table(comp, reliability) if reliability else None
    conf = None
    if rel is not None:
        conf = confidence_table(preds, reliability)
        paths["reliability"] = out_dir / f"{dataset_id}.reliability.tsv"
        paths["confidence"] = out_dir / f"{dataset_id}.confidence.tsv"
        rel.to_csv(paths["reliability"], sep="\t", index=False)
        conf.to_csv(paths["confidence"], sep="\t", index=False)
    consistency = {"status": "not_run", "reason": "no sample column found; pass --sample-key"}
    cons_table = None
    if sample_key and not cluster_key:
        # clusters computed here would follow the same unintegrated structure the refinement votes over, so calls
        # look consistent inside them by construction; only clusters from an integration of the samples can tell
        consistency = {"status": "not_assessable", "sample_key": sample_key,
                       "reason": "no usable cluster column; integrate the samples (e.g. Seurat or Harmony), "
                                 "cluster, and pass the cluster column with --cluster-key"}
    elif sample_key:
        samples = og_adata.obs[sample_key].to_numpy()
        clusters = og_adata.obs[cluster_key].to_numpy()
        match = cluster_correspondence(og_adata.X, samples, clusters)
        if match is not None and match < MIN_CLUSTER_CORRESPONDENCE:
            # the same cluster id stands for different cells in different samples, so comparing calls within a
            # cluster would report disagreement that comes from the clustering, not from the calls
            consistency = {"status": "not_assessable",
                           "reason": f"the clusters in {cluster_key!r} do not correspond across samples (only "
                                     f"{match:.0%} of clusters are most similar to the same cluster in the other "
                                     f"samples; integrated clusters of reference atlases: 99-100%); re-integrate "
                                     "the samples and recluster, or pass another column with --cluster-key"}
        else:
            cons_table, consistency = sample_consistency(preds["predicted_cell_type"].to_numpy(), samples, clusters)
            paths["sample_consistency"] = out_dir / f"{dataset_id}.sample_consistency.tsv"
            cons_table.to_csv(paths["sample_consistency"], sep="\t", index=False)
        consistency.update({"sample_key": sample_key, "cluster_key": cluster_key,
                            "cluster_correspondence": _round(match)})
    if sample_key:
        consistency["median_sample_composition_tv"] = _round(
            sample_composition_tv(preds["predicted_cell_type"].to_numpy(), og_adata.obs[sample_key].to_numpy()))
    write_html(paths["html"], f"scPCIP annotation report — {dataset_id}", meta or {}, comp, markers, figs,
               rel=rel, reliability=reliability, conf=conf, consistency=consistency, cons_table=cons_table)
    _dump_consistency(out_dir / f"{dataset_id}.sample_consistency.json", consistency)
    out = {k: str(v) for k, v in {**paths, **{f"fig_{k}": v for k, v in figs.items()}}.items()}
    out["sample_consistency_summary"] = consistency
    return out


def _round(x: float | None) -> float | None:
    return None if x is None else round(x, 4)


def _dump_consistency(path: Path, consistency: dict) -> None:
    path.write_text(json.dumps(consistency, indent=1), encoding="utf-8")


__all__ = ["top_variance_features", "embed", "composition", "marker_table", "group_scores", "make_figures",
           "make_extra_figures", "margin_quartile", "confidence_table", "reliability_table", "write_html",
           "build_report"]
