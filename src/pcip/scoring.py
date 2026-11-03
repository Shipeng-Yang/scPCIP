"""Rank-based program scoring over the cross-species orthogroup (OG) layer.

Scores use within-cell ranks, so they do not depend on library size (scaling a cell's counts leaves its
scores unchanged). No absolute count thresholds are used.

Four methods share one per-cell rank engine:
  * weighted_rank (default): weighted mean rank percentile of the positive features minus
    lambda times that of the negative features.
  * ucell_like: UCell-style signature score with a maxRank cap.
  * aucell_like: recovery-curve AUC of the positive set within the top of the cell ranking.
  * hybrid_rank_auc = 0.5*weighted_rank + 0.5*aucell_like.

Program OGs absent from the dataset (common across species) are not zero-filled, which would bias the
ranks; they are left out of the score and reported as coverage / n_missing_positive.

Output: a long ``program_scores`` table + ``obsm['X_pcip_scores']`` (cells x
programs primary score) with the program-id order stored in ``uns['pcip_score_programs']``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

SCORING_METHODS = ("weighted_rank", "ucell_like", "aucell_like", "hybrid_rank_auc")

_SCORE_COLUMNS = [
    "cell_id", "dataset_id", "species", "program_id", "method", "score",
    "positive_score", "negative_score", "coverage", "n_detected_positive",
    "n_missing_positive",
]


@dataclass
class ProgramFeatures:
    """A program's positive/negative features resolved against a dataset's var axis."""

    program_id: str
    pos_idx: np.ndarray          # var-column indices of positive OGs present in this dataset
    pos_w: np.ndarray            # matching weights (>=0), aligned to pos_idx
    neg_idx: np.ndarray
    neg_w: np.ndarray
    n_pos_total: int             # positive OGs in the dictionary (present + missing)
    n_neg_total: int

    @property
    def coverage(self) -> float:
        return self.pos_idx.size / self.n_pos_total if self.n_pos_total else 0.0

    @property
    def n_missing_positive(self) -> int:
        return self.n_pos_total - int(self.pos_idx.size)


def resolve_program_features(program: dict[str, Any], var_index: dict[str, int]) -> ProgramFeatures:
    """Map a program's OG feature ids to dataset var indices; absent OGs are dropped
    (counted missing, never zero-filled)."""
    def _resolve(feats: list[dict[str, Any]]):
        idx, w = [], []
        for f in feats or []:
            j = var_index.get(f.get("feature_id"))
            if j is not None:
                idx.append(j)
                w.append(float(f.get("weight", 1.0)))
        return np.asarray(idx, dtype=int), np.asarray(w, dtype=float)

    pos = program.get("positive_features", [])
    neg = program.get("negative_features", [])
    pos_idx, pos_w = _resolve(pos)
    neg_idx, neg_w = _resolve(neg)
    return ProgramFeatures(
        program_id=program.get("program_id", ""),
        pos_idx=pos_idx, pos_w=pos_w, neg_idx=neg_idx, neg_w=neg_w,
        n_pos_total=len(pos), n_neg_total=len(neg),
    )


def _percentile_ranks(block: np.ndarray) -> np.ndarray:
    """Per-row (cell) percentile rank in [0,1]; higher expression -> higher value.
    Ties get the average rank. Sparse zeros are real values tied at the bottom."""
    from scipy.stats import rankdata

    n = block.shape[1]
    out = np.empty_like(block, dtype=np.float64)
    for i in range(block.shape[0]):
        r = rankdata(block[i], method="average")     # 1..n, higher expr -> higher rank
        out[i] = (r - 1.0) / (n - 1) if n > 1 else np.zeros(n)
    return out


def _to_dense_block(X, rows: slice):
    import scipy.sparse as sp
    sub = X[rows]
    return sub.toarray() if sp.issparse(sub) else np.asarray(sub)


def _wmean(pct_block: np.ndarray, idx: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted mean of percentile ranks over feature columns idx (per cell)."""
    if idx.size == 0:
        return np.zeros(pct_block.shape[0])
    sub = pct_block[:, idx]                          # cells x |idx|
    sw = w.sum()
    return (sub * w).sum(axis=1) / sw if sw > 0 else sub.mean(axis=1)


def _desc_rank(pct_block: np.ndarray) -> np.ndarray:
    """Descending rank (1 = highest expression .. n = lowest) from ascending percentile."""
    n = pct_block.shape[1]
    return n - pct_block * (n - 1)            # pct=1 (top) -> rank 1; pct=0 -> rank n


def _ucell(pct_block: np.ndarray, idx: np.ndarray, maxrank_frac: float) -> np.ndarray:
    """Canonical UCell signature score (unweighted): descending ranks capped at a
    top maxRank window (genes below the window are tied at maxRank+1), normalized via
    the Mann-Whitney-U formula. maxRank cap is what distinguishes it from a
    plain weighted-rank mean."""
    if idx.size == 0:
        return np.zeros(pct_block.shape[0])
    n_feat = pct_block.shape[1]
    max_rank = max(1, int(round(maxrank_frac * n_feat)))
    desc = _desc_rank(pct_block)[:, idx]
    capped = np.minimum(desc, max_rank + 1)
    n_sig = idx.size
    sum_r = capped.sum(axis=1)
    u = sum_r - n_sig * (n_sig + 1) / 2.0
    score = 1.0 - u / (n_sig * max_rank)
    return np.clip(score, 0.0, 1.0)


def _aucell(pct_block: np.ndarray, idx: np.ndarray, top_fraction: float) -> np.ndarray:
    """AUCell-style recovery-curve AUC: area under the cumulative recovery of the
    positive set within the top `top_fraction` of the cell ranking, normalized to
    [0,1] (a positive gene higher in the ranking contributes more area)."""
    if idx.size == 0:
        return np.zeros(pct_block.shape[0])
    n_feat = pct_block.shape[1]
    thresh = max(1, int(round(top_fraction * n_feat)))
    desc = _desc_rank(pct_block)[:, idx]
    contrib = np.clip(thresh - desc + 1.0, 0.0, thresh)   # rank 1 -> thresh, rank>thresh -> 0
    return contrib.sum(axis=1) / (thresh * idx.size)


def _method_score(pct_block: np.ndarray, pf: ProgramFeatures, method: str,
                  lam: float, top_fraction: float, ucell_maxrank_frac: float
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (score, positive_score, negative_score) for one program on a rank block."""
    neg = _wmean(pct_block, pf.neg_idx, pf.neg_w)
    if method == "aucell_like":
        pos = _aucell(pct_block, pf.pos_idx, top_fraction)
    elif method == "ucell_like":
        pos = _ucell(pct_block, pf.pos_idx, ucell_maxrank_frac)
    elif method == "hybrid_rank_auc":
        pos = 0.5 * _wmean(pct_block, pf.pos_idx, pf.pos_w) + 0.5 * _aucell(pct_block, pf.pos_idx, top_fraction)
    else:  # weighted_rank (default)
        pos = _wmean(pct_block, pf.pos_idx, pf.pos_w)
    return pos - lam * neg, pos, neg


@dataclass
class ScoreResult:
    program_ids: list[str]
    scores: np.ndarray                       # cells x programs (primary score, method)
    method: str
    coverage: np.ndarray                     # per program (dataset-level)
    n_missing_positive: np.ndarray           # per program
    n_detected_positive: np.ndarray = field(repr=False, default=None)  # cells x programs
    positive_score: np.ndarray = field(repr=False, default=None)
    negative_score: np.ndarray = field(repr=False, default=None)


def score_cells(adata: Any, dictionary: dict[str, Any], method: str = "weighted_rank",
                lam: float = 1.0, top_fraction: float = 0.05, ucell_maxrank_frac: float = 0.05,
                chunk: int = 4000, write_obsm: bool = True) -> ScoreResult:
    """Score every cell against every program (library-size-invariant rank scores).

    Writes ``obsm['X_pcip_scores']`` (cells x programs) + ``uns['pcip_score_programs']``
    (id order) + ``uns['pcip_score_method']`` when write_obsm. Missing OGs are excluded
    and reported via coverage / n_missing_positive (never zero-filled)."""
    if method not in SCORING_METHODS:
        raise ValueError(f"unknown scoring method {method!r}; choose from {SCORING_METHODS}")
    if chunk <= 0:
        raise ValueError(f"chunk must be > 0, got {chunk}")
    if not (0.0 < top_fraction <= 1.0):
        raise ValueError(f"top_fraction must be in (0,1], got {top_fraction}")
    if not (0.0 < ucell_maxrank_frac <= 1.0):
        raise ValueError(f"ucell_maxrank_frac must be in (0,1], got {ucell_maxrank_frac}")
    programs = dictionary.get("programs", [])
    var_index = {str(v): i for i, v in enumerate(adata.var_names)}
    pfs = [resolve_program_features(p, var_index) for p in programs]
    n_cells, n_prog = adata.n_obs, len(programs)

    scores = np.zeros((n_cells, n_prog))
    pos_all = np.zeros((n_cells, n_prog))
    neg_all = np.zeros((n_cells, n_prog))
    ndet = np.zeros((n_cells, n_prog), dtype=int)

    X = adata.X
    for start in range(0, n_cells, chunk):
        rows = slice(start, min(start + chunk, n_cells))
        raw = _to_dense_block(X, rows)               # cells_block x n_var (raw counts)
        pct = _percentile_ranks(raw)
        for j, pf in enumerate(pfs):
            s, pos, neg = _method_score(pct, pf, method, lam, top_fraction, ucell_maxrank_frac)
            scores[rows, j] = s
            pos_all[rows, j] = pos
            neg_all[rows, j] = neg
            if pf.pos_idx.size:
                ndet[rows, j] = (raw[:, pf.pos_idx] > 0).sum(axis=1)

    coverage = np.array([pf.coverage for pf in pfs])
    n_missing = np.array([pf.n_missing_positive for pf in pfs])

    if write_obsm:
        adata.obsm["X_pcip_scores"] = scores
        adata.uns["pcip_score_programs"] = [pf.program_id for pf in pfs]
        adata.uns["pcip_score_method"] = method

    return ScoreResult(
        program_ids=[pf.program_id for pf in pfs], scores=scores, method=method,
        coverage=coverage, n_missing_positive=n_missing, n_detected_positive=ndet,
        positive_score=pos_all, negative_score=neg_all,
    )


def scores_to_long_df(adata: Any, result: ScoreResult, dataset_id: str = "", species: str = ""):
    """Flatten a ScoreResult to the long program_scores table."""
    import pandas as pd

    cell_ids = list(map(str, adata.obs_names))
    rows = []
    for j, pid in enumerate(result.program_ids):
        cov = float(result.coverage[j])
        nmiss = int(result.n_missing_positive[j])
        for i, cid in enumerate(cell_ids):
            rows.append((cid, dataset_id, species, pid, result.method,
                         round(float(result.scores[i, j]), 6),
                         round(float(result.positive_score[i, j]), 6),
                         round(float(result.negative_score[i, j]), 6),
                         round(cov, 4), int(result.n_detected_positive[i, j]), nmiss))
    return pd.DataFrame(rows, columns=_SCORE_COLUMNS)


def weighted_rank_score(adata: Any, program: dict[str, Any], lam: float = 1.0) -> np.ndarray:
    """Convenience: weighted_rank score of one program for all cells."""
    res = score_cells(adata, {"programs": [program]}, method="weighted_rank",
                       lam=lam, write_obsm=False)
    return res.scores[:, 0]


__all__ = [
    "SCORING_METHODS", "ProgramFeatures", "resolve_program_features",
    "ScoreResult", "score_cells", "scores_to_long_df", "weighted_rank_score",
]
