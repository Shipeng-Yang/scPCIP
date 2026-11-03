"""Multi-dataset program learning: cross-dataset consensus candidates + sparse supervised weights
with stability selection.

Every feature is the within-cell rank percentile of an orthogroup (the same representation the rank-based
scorer uses), so learned programs transfer to scaled deposits and other library sizes.

Algorithm (per coarse type t, one-vs-rest):
  1. candidates: per training dataset, contrast_d(f) = mean pct of f in t - mean over the other types' means;
                   keep OGs with contrast_d > 0 in at least ``min_datasets`` datasets that contain t, ranked by the
                   median contrast across those datasets; top ``n_candidates``.
  2. stability: ``n_rounds`` rounds; each round drops one random training dataset (when >= 3 contain t) and
                   draws a per-type cell bootstrap; fits L1 logistic regression on the candidates; records which
                   OGs get a positive coefficient and their size.
  3. program: OGs selected in >= ``min_freq`` of rounds, ordered by frequency then median coefficient;
                   at most ``max_features``; weight = median positive coefficient / its max (so weights lie in (0, 1]).
                   Fewer than ``min_features`` stable OGs -> status "few_stable" (kept as is, never padded);
                   none -> status "no_stable_features" and the type gets no program.
The output is a standard PCIP dictionary (positive_features with weights) usable by ``pcip.predict``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np


@dataclass
class LearnConfig:
    n_candidates: int = 300
    min_datasets: int = 2
    n_rounds: int = 30
    cells_per_type: int = 400
    l1_c: float = 0.2
    min_freq: float = 0.6
    max_features: int = 60
    min_features: int = 10
    seed: int = 20261003
    negatives: bool = False      # also stability-select negative features (genes OFF in the type)
    max_negative: int = 30
    balance_groups: bool = False  # per-type bootstrap split evenly across groups (species), so one
    #                               dominant group (e.g. one species with most cells) cannot set the coefficients alone
    extra: dict = field(default_factory=dict)


def rank_pct(X) -> np.ndarray:
    """Within-cell average-rank percentile (cells x features), as used by the rank-based scorer."""
    from scipy.stats import rankdata
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return rankdata(X, axis=1, method="average") / X.shape[1]


def type_contrast(pct: np.ndarray, labels: np.ndarray, target: str, types: Sequence[str]) -> np.ndarray:
    others = [t for t in types if t != target and np.any(labels == t)]
    ref = np.mean([pct[labels == r].mean(axis=0) for r in others], axis=0)
    return pct[labels == target].mean(axis=0) - ref


def consensus_candidates(blocks: Sequence[tuple[np.ndarray, np.ndarray]], target: str, types: Sequence[str],
                         n: int, min_datasets: int) -> np.ndarray:
    """Indices of candidate features for `target` from per-dataset (pct, labels) blocks."""
    # a block needs the target AND at least one comparator type; a target-only block has no contrast (NaN) and
    # would turn every feature's median NaN
    cons = [type_contrast(p, l, target, types) for p, l in blocks
            if np.any(l == target) and any(np.any(l == t) for t in types if t != target)]
    if not cons:
        return np.array([], dtype=int)
    C = np.vstack(cons)
    need = min(min_datasets, C.shape[0])
    ok = (C > 0).sum(axis=0) >= need
    med = np.median(C, axis=0)
    order = np.argsort(-med, kind="stable")
    return np.array([i for i in order if ok[i] and med[i] > 0][:n], dtype=int)


def negative_candidates(blocks: Sequence[tuple[np.ndarray, np.ndarray]], target: str, types: Sequence[str],
                        n: int, min_datasets: int) -> np.ndarray:
    """Indices of features consistently lower in `target` than in the other types (mirror of consensus_candidates)."""
    cons = [-type_contrast(p, l, target, types) for p, l in blocks
            if np.any(l == target) and any(np.any(l == t) for t in types if t != target)]
    if not cons:
        return np.array([], dtype=int)
    C = np.vstack(cons)
    need = min(min_datasets, C.shape[0])
    ok = (C > 0).sum(axis=0) >= need
    med = np.median(C, axis=0)
    order = np.argsort(-med, kind="stable")
    return np.array([i for i in order if ok[i] and med[i] > 0][:n], dtype=int)


def _bootstrap(labels: np.ndarray, groups: np.ndarray, keep_groups: set, cap: int, rng,
               balance: bool = False) -> np.ndarray:
    idx = []
    for t in np.unique(labels):
        pool = np.flatnonzero((labels == t) & np.isin(groups, list(keep_groups)))
        if not pool.size:
            continue
        if not balance:
            idx.append(rng.choice(pool, size=min(cap, pool.size), replace=True))
            continue
        gs = np.unique(groups[pool])            # each group holding type t gets an equal share of the cap
        per = max(1, min(cap, pool.size) // gs.size)
        for g in gs:
            sub = pool[groups[pool] == g]
            idx.append(rng.choice(sub, size=per, replace=True))
    return np.concatenate(idx)


def learn_programs(pct: np.ndarray, labels: np.ndarray, groups: np.ndarray, feats: Sequence[str],
                   types: Sequence[str], cfg: LearnConfig | None = None) -> tuple[dict, dict]:
    """Learn one positive program per type. `groups` = training dataset of each cell.
    Returns (dictionary, diagnostics)."""
    from sklearn.linear_model import LogisticRegression
    cfg = cfg or LearnConfig()
    rng = np.random.default_rng(cfg.seed)
    feats = np.asarray(list(map(str, feats)))
    gs = list(dict.fromkeys(groups.tolist()))
    blocks = [(pct[groups == g], labels[groups == g]) for g in gs]
    programs, diag = [], {}
    for t in types:
        if not np.any(labels == t):
            continue
        cand = consensus_candidates(blocks, t, types, cfg.n_candidates, cfg.min_datasets)
        if cand.size < cfg.min_features:
            diag[t] = {"status": "too_few_candidates", "n_candidates": int(cand.size)}
            continue
        n_pos = cand.size
        if cfg.negatives:                 # fit on positive + negative candidates; positives are still read from cand
            neg = negative_candidates(blocks, t, types, cfg.n_candidates, cfg.min_datasets)
            cand = np.concatenate([cand, neg[~np.isin(neg, cand)]])
        has_t = [g for g in gs if np.any((groups == g) & (labels == t))]
        hits = np.zeros(cand.size)
        coefs = [[] for _ in range(cand.size)]
        nhits = np.zeros(cand.size)
        ncoefs = [[] for _ in range(cand.size)]
        for _ in range(cfg.n_rounds):
            keep = set(gs)
            if len(has_t) >= 3:
                keep.discard(has_t[int(rng.integers(len(has_t)))])
            idx = _bootstrap(labels, groups, keep, cfg.cells_per_type, rng, balance=cfg.balance_groups)
            y = (labels[idx] == t).astype(int)
            if y.min() == y.max():
                continue
            m = LogisticRegression(penalty="l1", C=cfg.l1_c, solver="liblinear", class_weight="balanced", random_state=cfg.seed,
                                   max_iter=2000)
            m.fit(pct[np.ix_(idx, cand)], y)
            w = m.coef_[0]
            for j in np.flatnonzero(w > 0):
                hits[j] += 1
                coefs[j].append(w[j])
            for j in np.flatnonzero(w < 0):
                nhits[j] += 1
                ncoefs[j].append(-w[j])
        freq = hits / cfg.n_rounds
        med = np.array([np.median(c) if c else 0.0 for c in coefs])
        sel = [j for j in np.lexsort((-med, -freq)) if j < n_pos and freq[j] >= cfg.min_freq and med[j] > 0
               ][:cfg.max_features]
        if not sel:                       # nothing passed stability selection: report, never fill with unstable OGs
            diag[t] = {"status": "no_stable_features", "n_candidates": int(cand.size), "n_selected": 0,
                       "n_freq_ge_min": int((freq >= cfg.min_freq).sum())}
            continue
        wmax = max(med[j] for j in sel)
        prog = {"program_id": f"PCIP_ATLAS_ROOT_{t.upper()}", "cell_type_coarse": t,
                "program_type": "identity", "scope": "root", "level": "coarse",
                "positive_features": [{"feature_id": str(feats[cand[j]]),
                                       "weight": round(float(med[j] / wmax), 4)} for j in sel]}
        if cfg.negatives:
            nfreq = nhits / cfg.n_rounds
            nmed = np.array([np.median(c) if c else 0.0 for c in ncoefs])
            nsel = [j for j in np.lexsort((-nmed, -nfreq)) if nfreq[j] >= cfg.min_freq and nmed[j] > 0
                    ][:cfg.max_negative]
            if nsel:
                nmax = max(nmed[j] for j in nsel)
                prog["negative_features"] = [{"feature_id": str(feats[cand[j]]),
                                              "weight": round(float(nmed[j] / nmax), 4)} for j in nsel]
        programs.append(prog)
        diag[t] = {"status": "ok" if len(sel) >= cfg.min_features else "few_stable",
                   "n_candidates": int(cand.size), "n_selected": len(sel),
                   "min_freq_selected": float(min(freq[j] for j in sel)),
                   "n_freq_ge_min": int((freq >= cfg.min_freq).sum())}
    return {"version": "atlas_root_learned", "feature_layer": "orthogroup", "programs": programs}, diag


def consensus_programs(pct: np.ndarray, labels: np.ndarray, groups: np.ndarray, feats: Sequence[str],
                       types: Sequence[str], n: int = 45, min_datasets: int = 2) -> dict:
    """Ablation arm: top-n consensus-contrast OGs per type, uniform weights (no supervised model)."""
    feats = np.asarray(list(map(str, feats)))
    gs = list(dict.fromkeys(groups.tolist()))
    blocks = [(pct[groups == g], labels[groups == g]) for g in gs]
    progs = []
    for t in types:
        if not np.any(labels == t):
            continue
        cand = consensus_candidates(blocks, t, types, n, min_datasets)
        if cand.size:
            progs.append({"program_id": f"PCIP_CONSENSUS_ROOT_{t.upper()}", "cell_type_coarse": t,
                          "program_type": "identity",
                          "positive_features": [{"feature_id": str(feats[i]), "weight": 1.0} for i in cand]})
    return {"version": "atlas_root_consensus", "feature_layer": "orthogroup", "programs": progs}


__all__ = ["LearnConfig", "rank_pct", "type_contrast", "consensus_candidates", "learn_programs",
           "consensus_programs"]
