"""Hierarchical program dictionaries: coarse cell type, then fine type within it, with fallback.

A hierarchical dictionary is a standard PCIP dictionary (``programs`` = coarse programs, usable by
``pcip.predict`` as is) plus ``fine_programs``: {coarse type: [fine programs]}. Fine programs are learned only
from cells that carry a fine label of that coarse type (cells labelled at coarse level only train the coarse
programs), and only when at least two fine types of that coarse type have labelled cells.

Prediction: argmax over coarse programs; then, if the predicted coarse type has fine programs, argmax over those.
The fine call is kept only when its margin (top minus runner-up fine score) is >= ``fine_margin``; otherwise the cell
is reported at coarse level. Output columns: predicted_coarse, coarse_margin, predicted_fine ('' when none),
fine_margin (NaN when none), predicted_label (fine if kept, else coarse), predicted_level ('fine' / 'coarse').

Score correction (coarse level, uses only the query's own scores, never labels): ``score_correction="bg_center"``
subtracts from each program's scores its median over the query cells NOT called that type by the raw argmax (the
program's background in this query); a program that wins no cell is centred on its median over all cells. This
removes "sink" programs whose baseline is high on every cell of a distant species (without it, 75-96% of monocot
root cells were called meristem / exodermis). "center" subtracts the all-cell median; "colz" z-scores each program
across cells; "none" (default) leaves scores unchanged.

Ubiquity filter (``ubiquity_max``, label-free, off by default): features detected (count > 0) in at least that
fraction of the query's cells are dropped from every program's positive features before scoring (a program that
would lose all of them is kept whole). Example: orthogroups only mildly pollen-biased in Brassicaceae buds sit in
the top 10-30% of almost every rice nucleus and made the pollen program win 26-64% of rice cells.

Query refinement (``refine="selftrain_vote"``, also label-free): (1) per predicted coarse type, the cells in the top
quartile of coarse margin become pseudo-labels; a multinomial logistic regression on the query's own PCA embedding
(normalise / log1p unless scaled, top-2000-variance OGs, 30 PCs) is fitted on them and relabels every cell;
(2) Leiden clusters (resolution 2.0) whose majority label covers >= 50% are set to that label. Fine types are then
called inside the final coarse type. Leave-one-species-out over 12 species: coarse macro-F1 species mean 0.336
(bg_center) -> 0.406.
"""
from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from . import program_learn as pl


def learn_hierarchy(pct: np.ndarray, coarse: np.ndarray, fine: np.ndarray, groups: np.ndarray,
                    feats: Sequence[str], cfg: pl.LearnConfig | None = None,
                    state_types: Sequence[str] = ()) -> tuple[dict, dict]:
    """`coarse` / `fine` = per-cell labels ('' = unlabelled at that level). Returns (dictionary, diagnostics).
    `state_types` (e.g. meristem): coarse labels that are developmental states, not lineages, learned in
    the coarse step together with the lineages (their cells are contrast for the lineage programs), then stored
    under ``state_programs``: never predicted as a cell type, reported as a state score."""
    coarse, fine, groups = (np.asarray(v).astype(str) for v in (coarse, fine, groups))
    lab = coarse != ""
    ctypes = sorted(set(coarse[lab]))
    d, cdiag = pl.learn_programs(pct[lab], coarse[lab], groups[lab], feats, ctypes, cfg)
    fine_programs, fdiag = {}, {}
    for c in ctypes:
        if c in state_types:
            continue
        m = (coarse == c) & (fine != "")
        ftypes = sorted(set(fine[m]))
        if len(ftypes) < 2:
            continue
        fd, fdi = pl.learn_programs(pct[m], fine[m], groups[m], feats, ftypes, cfg)
        progs = []
        for p in fd["programs"]:
            f = p["cell_type_coarse"]
            progs.append({**p, "program_id": f"PCIP_ATLAS_ROOT_{c.upper()}__{f.upper()}".replace(" ", "_"),
                          "cell_type_coarse": c, "cell_type_fine": f, "parent_coarse": c, "level": "fine"})
        fdiag[c] = {"n_cells": int(m.sum()), "fine_types": ftypes, "per_type": fdi}
        if len(progs) >= 2:
            fine_programs[c] = progs
    d = {**d, "version": "atlas_root_hierarchical", "fine_programs": fine_programs}
    if state_types:
        d["state_programs"] = [p for p in d["programs"] if p["cell_type_coarse"] in state_types]
        d["programs"] = [p for p in d["programs"] if p["cell_type_coarse"] not in state_types]
    return d, {"coarse": cdiag, "fine": fdiag}


SCORE_CORRECTIONS = ("none", "bg_center", "center", "colz")


def correct_scores(M: np.ndarray, method: str = "none") -> np.ndarray:
    """Per-query correction of a cells x types score matrix (see module docstring)."""
    if method not in SCORE_CORRECTIONS:
        raise ValueError(f"score_correction must be one of {SCORE_CORRECTIONS}, got {method!r}")
    M = np.asarray(M, dtype=float)
    if method == "none":
        return M
    if method == "center":
        return M - np.nanmedian(M, axis=0)
    if method == "colz":
        return (M - np.nanmean(M, axis=0)) / (np.nanstd(M, axis=0) + 1e-9)
    raw = np.nanargmax(np.where(np.isnan(M), -np.inf, M), axis=1)
    out = M.copy()
    for j in range(M.shape[1]):
        others = raw != j
        out[:, j] = M[:, j] - np.nanmedian(M[others, j] if others.any() else M[:, j])
    return out


def cluster_vote(labels: np.ndarray, clusters: np.ndarray, min_purity: float = 0.0) -> np.ndarray:
    """Replace each cell's label by the unique majority label of its cluster (a tie at the top is not voted).

    Only clusters whose majority label covers >= ``min_purity`` of the cluster are voted; mixed clusters keep the
    per-cell labels (guards against a coarse cluster that merged several cell types being relabelled as one)."""
    labels, clusters = np.asarray(labels).astype(str), np.asarray(clusters).astype(str)
    out = labels.copy()
    for c in np.unique(clusters):
        m = clusters == c
        vals, counts = np.unique(labels[m], return_counts=True)
        unique_top = (counts == counts.max()).sum() == 1          # a tie at the top has no majority to vote for
        if unique_top and counts.max() / m.sum() >= min_purity:
            out[m] = vals[np.argmax(counts)]
    return out


def cap_refine_growth(unrefined: np.ndarray, refined: np.ndarray, cap: float | None) -> np.ndarray:
    """Refine growth cap: a type that refinement grows beyond ``cap`` x its pre-refinement cell count
    gets none of the cells moved into it (they revert to their pre-refinement label); cells leaving a type are
    untouched. Guards against refinement inventing a type the query lacks (an shr-2 mutant root: endodermis 1.9%
    of cells before refinement -> 10.9% after). ``cap=None`` returns ``refined`` unchanged."""
    unrefined = np.asarray(unrefined).astype(str).astype(object)    # object: reverted labels may be longer than
    refined = np.asarray(refined).astype(str).astype(object)        # any refined one (fixed-width <U would truncate)
    if cap is None:
        return refined.copy()
    out = refined.copy()
    for t in np.unique(refined):
        if (refined == t).sum() > cap * (unrefined == t).sum():
            moved_in = (refined == t) & (unrefined != t)
            out[moved_in] = unrefined[moved_in]
    return out


REFINES = ("none", "selftrain_vote")


def self_train(pca: np.ndarray, labels: np.ndarray, margin: np.ndarray, q: float = 0.75, seed: int = 0,
               min_cells: int = 10) -> np.ndarray:
    """Relabel cells with a classifier fitted on each label's top-(1-q) margin cells (transductive, label-free)."""
    from sklearn.linear_model import LogisticRegression
    labels = np.asarray(labels).astype(str)
    margin = np.asarray(margin, dtype=float)
    keep = np.zeros(len(labels), dtype=bool)
    for t in np.unique(labels):
        m = (labels == t) & np.isfinite(margin)
        if m.sum() < min_cells:
            continue
        keep |= m & (margin >= np.quantile(margin[m], q))
    if len(np.unique(labels[keep])) < 2:
        return labels.copy()
    clf = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced", random_state=seed)
    clf.fit(pca[keep], labels[keep])
    return clf.predict(pca).astype(str)


def query_embedding(og_adata: Any, resolution: float = 2.0, n_pcs: int = 30, seed: int = 0,
                    regress_depth: bool = False):
    """(PCA coordinates, Leiden cluster labels) of the query itself. ``regress_depth``: before PCA,
    regress log(1 + detected features per cell) out of every selected feature, so clusters cannot split by
    library complexity (otherwise refinement can turn depth quintiles into cell types)."""
    import scanpy as sc
    from . import omg_baseline as omg
    a = og_adata.copy()
    Xraw = a.X
    n_det = np.asarray((Xraw != 0).sum(axis=1)).ravel() if hasattr(Xraw, "toarray") else (np.asarray(Xraw) != 0).sum(1)
    sc.pp.filter_genes(a, min_cells=3)
    omg._normalise(a)
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    hv = np.argsort(-X.var(axis=0), kind="stable")[:min(2000, a.n_vars)]
    a = a[:, hv].copy()
    if regress_depth:
        Y = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X, dtype=float)
        d = np.log1p(n_det.astype(float))
        d -= d.mean()
        Y = Y - Y.mean(axis=0)
        a.X = Y - np.outer(d, (d @ Y) / max(float(d @ d), 1e-12))
    sc.pp.pca(a, n_comps=min(n_pcs, a.n_vars - 1, a.n_obs - 1), random_state=seed)
    sc.pp.neighbors(a, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return np.asarray(a.obsm["X_pca"]), a.obs["leiden"].astype(str).to_numpy()


def _zero_tie_pct(og_adata: Any) -> np.ndarray:
    """Per cell, the within-cell percentile rank a zero-count feature gets (ties at the average rank), i.e. the
    floor every unexpressed feature sits on in scoring._percentile_ranks."""
    import scipy.sparse as sp
    X = og_adata.X
    n = og_adata.n_vars
    nz = np.asarray((X != 0).sum(axis=1)).ravel() if sp.issparse(X) else (np.asarray(X) != 0).sum(axis=1)
    z = n - nz
    return np.where(z > 0, (z - 1) / 2.0 / max(n - 1, 1), 0.0)


def _negative_penalty(og_adata: Any, programs: list, lam: float = 0.5, **kw: Any):
    """Per program with ``negative_features``, the weighted-mean percentile of those features above the
    cell's zero floor (0 when none is expressed). The scorer's own ``pos - lam*neg`` would charge every
    cell the zero floor (~0.3-0.5) and add a library-complexity bias, so negatives are scored here instead.
    Returns {cell_type_coarse: penalty array} (empty when no program has negatives)."""
    from . import predict as pr
    neg = [{"program_id": p.get("program_id", ""), "cell_type_coarse": p["cell_type_coarse"],
            "positive_features": p["negative_features"]} for p in programs if p.get("negative_features")]
    if not neg:
        return {}
    # the floor is a percentile-rank quantity -> score negatives with the weighted-rank mean whatever the positive
    # scoring method is
    _, ns, _ = pr.predict_og_adata(og_adata, {"programs": neg}, **{**kw, "method": "weighted_rank"})
    floor = _zero_tie_pct(og_adata)
    return {c[len("score_"):]: lam * np.clip(ns[c].to_numpy(dtype=float) - floor, 0.0, None) for c in ns.columns}


def _strip_negatives(d: dict) -> dict:
    return {**d, "programs": [{k: v for k, v in p.items() if k != "negative_features"} for p in d.get("programs", [])]}


def drop_ubiquitous(og_adata: Any, hdict: dict, ubiquity_max: float | None) -> dict:
    """Copy of `hdict` without the positive features detected in >= `ubiquity_max` of the query cells (coarse, fine
    and state programs); a program that would keep no positive feature is left unchanged. None -> `hdict` itself."""
    if ubiquity_max is None:
        return hdict
    X = og_adata.X
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    ubiq = {str(v) for v, d in zip(og_adata.var_names, det) if d >= ubiquity_max}

    def prog(p: dict) -> dict:
        keep = [f for f in p.get("positive_features", []) if str(f["feature_id"]) not in ubiq]
        return {**p, "positive_features": keep} if keep else p
    out = {**hdict, "programs": [prog(p) for p in hdict.get("programs", [])]}
    if hdict.get("fine_programs"):
        out["fine_programs"] = {c: [prog(p) for p in ps] for c, ps in hdict["fine_programs"].items()}
    if hdict.get("state_programs"):
        out["state_programs"] = [prog(p) for p in hdict["state_programs"]]
    return out


def predict_hierarchy(og_adata: Any, hdict: dict, fine_margin: float = 0.0, score_correction: str = "none",
                      refine: str = "none", vote_resolution: float = 2.0, vote_min_purity: float = 0.5,
                      return_scores: bool = False, negative_lambda: float = 0.5, return_fine_scores: bool = False,
                      embed_regress_depth: bool = False, refine_growth_cap: float | None = None,
                      ubiquity_max: float | None = None, **kw: Any):
    """Hierarchical single prediction for an OG-layer AnnData; returns a DataFrame (one row per cell), plus the
    corrected coarse score matrix when ``return_scores``. ``coarse_margin`` always describes the final coarse label
    against those scores: score(final) - best other score (negative when refinement overrode the argmax);
    ``coarse_margin_unrefined`` keeps the pre-refinement argmax margin. ``return_fine_scores`` adds one column
    ``fine_score_<coarse>__<fine>`` per fine program (the score used for the fine call; NaN for cells not called
    that coarse type) and changes nothing else. ``ubiquity_max``: see the module docstring (None = off)."""
    import pandas as pd
    from . import predict as pr
    hdict = drop_ubiquitous(og_adata, hdict, ubiquity_max)
    coarse_dict = {k: v for k, v in hdict.items() if k not in ("fine_programs", "state_programs")}
    pen = _negative_penalty(og_adata, coarse_dict.get("programs", []), lam=negative_lambda, **kw)
    cp, cs, _ = pr.predict_og_adata(og_adata, _strip_negatives(coarse_dict), **kw)
    for t, v in pen.items():
        cs[f"score_{t}"] = cs[f"score_{t}"] - v
    ctypes = [c[len("score_"):] for c in cs.columns]
    M = correct_scores(cs.to_numpy(dtype=float), score_correction)
    if score_correction != "none" or pen:
        sp = pr.single_prediction(M, ctypes, seed=kw.get("seed", pr.DEFAULT_SEED))
        cp = cp.assign(predicted_cell_type=sp.predicted, margin=sp.margin, tie_at_top=sp.tie_at_top)
    if refine not in REFINES:
        raise ValueError(f"refine must be one of {REFINES}, got {refine!r}")
    out = pd.DataFrame({"cell_id": cp["cell_id"], "predicted_coarse": cp["predicted_cell_type"].astype(str),
                        "coarse_margin": cp["margin"].astype(float), "coarse_tie_at_top": cp["tie_at_top"].astype(bool)})
    if refine == "selftrain_vote":
        pca, clusters = query_embedding(og_adata, resolution=vote_resolution, seed=0, regress_depth=embed_regress_depth)
        st = self_train(pca, out["predicted_coarse"].to_numpy(), out["coarse_margin"].to_numpy(),
                        seed=kw.get("seed", 0) or 0)
        out["predicted_coarse_unrefined"] = out["predicted_coarse"]
        out["coarse_margin_unrefined"] = out["coarse_margin"]
        out["predicted_coarse"] = cap_refine_growth(out["predicted_coarse_unrefined"].to_numpy(),
                                                    cluster_vote(st, clusters, min_purity=vote_min_purity),
                                                    refine_growth_cap)
        # margin / tie of the final label against the corrected scores
        idx = {t: k for k, t in enumerate(ctypes)}
        j = np.array([idx.get(t, -1) for t in out["predicted_coarse"]])
        ok = j >= 0
        own = np.full(len(j), np.nan)
        own[ok] = M[np.flatnonzero(ok), j[ok]]
        masked = np.where(np.arange(M.shape[1])[None, :] == j[:, None], -np.inf, M)
        mg = own - np.max(masked, axis=1)
        out["coarse_margin"] = mg
        out["coarse_tie_at_top"] = np.isclose(mg, 0.0) & ok
    out["predicted_fine"] = ""
    out["fine_margin"] = np.nan
    fp = hdict.get("fine_programs") or {}
    if fp:
        flat = [{**p, "cell_type_coarse": f"{c}::{p['cell_type_fine']}"} for c, ps in fp.items() for p in ps]
        fpen = _negative_penalty(og_adata, flat, lam=negative_lambda, **kw)
        _, fs, _ = pr.predict_og_adata(og_adata, _strip_negatives({"programs": flat}), **kw)
        for t, v in fpen.items():
            fs[f"score_{t}"] = fs[f"score_{t}"] - v
        cols = [c[len("score_"):] for c in fs.columns]
        S = fs.to_numpy(dtype=float)
        if return_fine_scores:              # fixed schema: one column per fine program even if no cell is called
            for n in cols:                  # its coarse type
                cc, ff = n.split("::", 1)
                out[f"fine_score_{cc}__{ff}"] = np.nan
        for c in fp:
            j = [k for k, n in enumerate(cols) if n.split("::", 1)[0] == c]
            rows = np.flatnonzero(out["predicted_coarse"].to_numpy() == c)
            if not rows.size or len(j) < 2:
                continue
            sub = S[np.ix_(rows, j)]
            order = np.argsort(-sub, axis=1, kind="stable")
            top = sub[np.arange(rows.size), order[:, 0]]
            ru = sub[np.arange(rows.size), order[:, 1]]
            names = np.array([cols[k].split("::", 1)[1] for k in j])
            out.loc[out.index[rows], "predicted_fine"] = names[order[:, 0]]
            out.loc[out.index[rows], "fine_margin"] = top - ru
            if return_fine_scores:
                for jj, k in enumerate(j):
                    out.loc[out.index[rows], f"fine_score_{c}__{names[jj]}"] = sub[:, jj]
    sp_ = hdict.get("state_programs") or []
    if sp_:
        # state score = the state program's raw score, a ranking within the query
        spen = _negative_penalty(og_adata, sp_, lam=negative_lambda, **kw)
        _, ss, _ = pr.predict_og_adata(og_adata, _strip_negatives({"programs": sp_}), **kw)
        for c in ss.columns:
            t = c[len("score_"):]
            out[f"state_{t}_score"] = ss[c].to_numpy(dtype=float) - spen.get(t, 0.0)
    keep = (out["predicted_fine"] != "") & (out["fine_margin"].fillna(-np.inf) >= fine_margin)
    out["predicted_label"] = np.where(keep, out["predicted_fine"], out["predicted_coarse"])
    out["predicted_level"] = np.where(keep, "fine", "coarse")
    if return_scores:
        return out, pd.DataFrame(M, index=cs.index, columns=cs.columns)
    return out
