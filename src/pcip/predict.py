"""Single-prediction cell-type annotation (program-first).

One predicted cell type per cell, the argmax of gene-program scores, with confidence kept as extra columns
(top score, runner-up type and score, margin, tie flag). Nothing here abstains; the confidence columns let the
user filter. Scoring uses `pcip.scoring.score_cells`.
When several programs share a cell type, the type score is the max over its programs.
Ties at the top are broken uniformly at random with a fixed seed and flagged (`tie_at_top`).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

DEFAULT_SEED = 20261002


def program_cell_types(dictionary: dict[str, Any]) -> list[str]:
    """Coarse cell type of each program, in dictionary order."""
    out = []
    for p in dictionary.get("programs", []):
        ct = p.get("cell_type_coarse") or p.get("cell_type") or p.get("cell_type_fine") or p.get("program_id")
        out.append(str(ct))
    return out


def collapse_to_types(scores: np.ndarray, program_types: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    """(cells x programs) -> (cells x unique types) by taking the max over programs of a type."""
    scores = np.asarray(scores, dtype=float)
    types = list(dict.fromkeys(program_types))
    idx = {t: [j for j, pt in enumerate(program_types) if pt == t] for t in types}
    out = np.column_stack([scores[:, idx[t]].max(axis=1) for t in types]) if types else np.zeros((scores.shape[0], 0))
    return out, types


@dataclass
class SinglePrediction:
    predicted: np.ndarray        # str, one label per cell
    top_score: np.ndarray
    runner_up: np.ndarray        # str ("" if only one type)
    runner_up_score: np.ndarray
    margin: np.ndarray           # top_score - runner_up_score (0 when tied)
    tie_at_top: np.ndarray       # bool
    invalid: np.ndarray          # bool: no finite score for this cell -> predicted ""


def single_prediction(type_scores: np.ndarray, types: Sequence[str], seed: int = DEFAULT_SEED,
                      atol: float = 1e-12) -> SinglePrediction:
    """Argmax over types; ties at the top broken uniformly at random (fixed seed) and flagged."""
    s = np.asarray(type_scores, dtype=float)
    n, k = s.shape
    if k == 0:
        raise ValueError("no cell types to predict")
    types = np.asarray(list(types), dtype=object)
    rng = np.random.default_rng(seed)
    invalid = ~np.isfinite(s).any(axis=1)
    s = np.where(np.isfinite(s), s, -np.inf)
    mx = s.max(axis=1, keepdims=True)
    at_top = np.isclose(s, mx, rtol=0, atol=atol)
    n_top = at_top.sum(axis=1)
    choice = np.argmax(s, axis=1)
    for i in np.flatnonzero(n_top > 1):
        choice[i] = rng.choice(np.flatnonzero(at_top[i]))
    top = s[np.arange(n), choice]
    if k > 1:
        rest = s.copy()
        rest[np.arange(n), choice] = -np.inf
        ru = np.argmax(rest, axis=1)
        ru_score = rest[np.arange(n), ru]
        ru_label = types[ru]
    else:
        ru_score = np.full(n, np.nan)
        ru_label = np.full(n, "", dtype=object)
    ru_score = np.where(np.isfinite(ru_score), ru_score, np.nan)
    predicted = types[choice].astype(object)
    predicted[invalid] = ""
    top = np.where(invalid, np.nan, top)
    margin = np.where(invalid, np.nan, top - np.nan_to_num(ru_score, nan=top))
    return SinglePrediction(predicted=predicted, top_score=top, runner_up=ru_label,
                            runner_up_score=ru_score, margin=margin,
                            tie_at_top=(n_top > 1) & ~invalid, invalid=invalid)


def predict_og_adata(og_adata: Any, dictionary: dict[str, Any], method: str = "weighted_rank",
                     seed: int = DEFAULT_SEED, **scoring_params: Any):
    """Score an OG-layer AnnData and return (predictions DataFrame, per-type score DataFrame)."""
    import pandas as pd
    from . import scoring
    res = scoring.score_cells(og_adata, dictionary, method=method, write_obsm=False, **scoring_params)
    type_scores, types = collapse_to_types(res.scores, program_cell_types(dictionary))
    sp = single_prediction(type_scores, types, seed=seed)
    cells = list(map(str, og_adata.obs_names))
    pred = pd.DataFrame({"cell_id": cells, "predicted_cell_type": sp.predicted, "top_score": sp.top_score,
                         "runner_up_cell_type": sp.runner_up, "runner_up_score": sp.runner_up_score,
                         "margin": sp.margin, "tie_at_top": sp.tie_at_top, "invalid_score": sp.invalid})
    scores = pd.DataFrame(type_scores, index=cells, columns=[f"score_{t}" for t in types])
    return pred, scores, res
