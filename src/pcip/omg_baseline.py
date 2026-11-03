"""OMG-style cluster-level marker overlap (used by the query refinement and for comparison).

Follows Chau et al., Nat Commun 2025 (PMC11695714), Methods:
  1. cluster the query (Louvain/Leiden on PCA), top-N (N=200) positive markers per cluster
     (Wilcoxon; min.pct 0.3, logFC >= 0.45; ranked by log fold-change);
  2. map markers to orthogroups (here: inputs are already OG-level matrices, so markers ARE OGs);
  3. for every (query cluster I, reference type J): X = shared OGs, Y = sum_j!=J OMG[I,j],
     Z = sum_i!=I OMG[i,J], W = total - X - Y - Z; Fisher exact test, upper tail; BH-FDR over all
     pairs; significant at FDR < 0.01;
  4. cluster label = reference type with the smallest FDR; every cell inherits its cluster's label.
Deviations (documented, for a single-prediction comparison): a cluster with no significant
reference still gets its best (smallest-FDR) type and is flagged `omg_significant=False`;
markers are computed on OG-aggregated matrices; scanpy replaces Seurat; scaled matrices (negative
values, e.g. SCT / scale.data deposits) are used as-is instead of normalise + log1p, and their markers are
ranked by the Wilcoxon statistic (log fold-change / min.pct are undefined on centred values).
"""
from __future__ import annotations

from typing import Any, Sequence

import numpy as np

N_MARKERS, MIN_PCT, MIN_LOGFC, FDR_ALPHA = 200, 0.3, 0.45, 0.01


def is_scaled(adata: Any) -> bool:
    """True when the matrix holds scaled / centred values (negatives present) rather than counts."""
    X = adata.X
    m = X.min() if hasattr(X, "min") else np.min(X)
    return bool(np.asarray(m).item() < 0)


def _normalise(a: Any) -> None:
    """In place: library-size normalise + log1p for count matrices; scaled matrices are used as-is."""
    import scanpy as sc
    if not is_scaled(a):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)


def top_markers(adata: Any, groupby: str, n: int = N_MARKERS, min_pct: float = MIN_PCT,
                min_logfc: float = MIN_LOGFC) -> dict[str, list[str]]:
    """Top-n positive markers per group (Wilcoxon one-vs-rest), filtered as in OMG, ranked by logFC."""
    import scanpy as sc
    a = adata.copy()
    scaled = is_scaled(a)
    _normalise(a)
    sc.tl.rank_genes_groups(a, groupby, method="wilcoxon", pts=True, use_raw=False)
    out = {}
    for g in a.obs[groupby].astype(str).unique():
        df = sc.get.rank_genes_groups_df(a, group=g)
        if scaled:
            # log fold-change and fraction-expressing are undefined on centred values: rank by the
            # Wilcoxon statistic and keep positive (up-in-group) features only
            df = df[df["scores"].to_numpy() > 0].sort_values("scores", ascending=False)
        else:
            pts = a.uns["rank_genes_groups"]["pts"][g].reindex(df["names"]).to_numpy()
            keep = (df["logfoldchanges"].to_numpy() >= min_logfc) & (pts >= min_pct)
            df = df[keep].sort_values("logfoldchanges", ascending=False)
        out[str(g)] = list(map(str, df["names"].head(n)))
    return out


def overlap_matrix(query: dict[str, Sequence[str]], ref: dict[str, Sequence[str]]):
    """(query groups x reference groups) counts of shared markers."""
    qk, rk = list(query), list(ref)
    rs = {r: set(ref[r]) for r in rk}
    m = np.array([[len(set(query[q]) & rs[r]) for r in rk] for q in qk], dtype=int)
    return m, qk, rk


def omg_fisher(overlap: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Upper-tail Fisher p-values with the OMG contingency, and BH-FDR over all pairs."""
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    o = np.asarray(overlap, dtype=int)
    total, row, col = o.sum(), o.sum(axis=1), o.sum(axis=0)
    p = np.ones_like(o, dtype=float)
    for i in range(o.shape[0]):
        for j in range(o.shape[1]):
            x = o[i, j]
            y, z = row[i] - x, col[j] - x
            w = total - x - y - z
            p[i, j] = fisher_exact([[x, y], [z, w]], alternative="greater")[1]
    fdr = multipletests(p.ravel(), method="fdr_bh")[1].reshape(p.shape)
    return p, fdr


def assign_clusters(fdr: np.ndarray, ref_types: Sequence[str], alpha: float = FDR_ALPHA):
    """Best (smallest-FDR) reference type per query cluster + significance flag."""
    best = np.argmin(fdr, axis=1)
    return [str(ref_types[b]) for b in best], [bool(fdr[i, b] < alpha) for i, b in enumerate(best)]


def cluster_query(adata: Any, resolution: float = 0.5, n_pcs: int = 30, seed: int = 0) -> np.ndarray:
    """Standard scanpy clustering (normalize, log1p, HVG, PCA, kNN, Leiden)."""
    import scanpy as sc
    a = adata.copy()
    sc.pp.filter_genes(a, min_cells=3)       # all-zero OGs (small gene panels) are uninformative
    _normalise(a)
    # top-variance features on log data (binned HVG flavours fail on skewed small-panel OG matrices)
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    var = X.var(axis=0)
    hv = np.zeros(a.n_vars, dtype=bool)
    hv[np.argsort(-var, kind="stable")[:min(2000, a.n_vars)]] = True
    a.var["highly_variable"] = hv
    sc.pp.pca(a, n_comps=min(n_pcs, int(hv.sum()) - 1, a.n_obs - 1), mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(a, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return a.obs["leiden"].astype(str).to_numpy()


def omg_predict(query_adata: Any, ref_markers: dict[str, Sequence[str]], resolution: float = 0.5,
                seed: int = 0, clusters: np.ndarray | None = None):
    """Per-cell labels from OMG cluster assignment. Returns (labels, clusters, info dict)."""
    q = query_adata.copy()
    q.obs["_omg_cluster"] = clusters if clusters is not None else cluster_query(q, resolution, seed=seed)
    qm = top_markers(q, "_omg_cluster")
    ov, qk, rk = overlap_matrix(qm, ref_markers)
    _, fdr = omg_fisher(ov)
    lab, sig = assign_clusters(fdr, rk)
    lab_of, sig_of = dict(zip(qk, lab)), dict(zip(qk, sig))
    cl = q.obs["_omg_cluster"].astype(str).to_numpy()
    labels = np.array([lab_of[c] for c in cl], dtype=object)
    info = {"n_clusters": len(qk), "cluster_labels": lab_of, "cluster_significant": sig_of,
            "overlap": ov.tolist(), "query_clusters": qk, "reference_types": rk}
    return labels, cl, info
