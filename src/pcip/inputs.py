"""Input checks and the gene -> orthogroup map builder.

`ensure_counts` makes sure scoring sees raw counts: many h5ad files exported from Seurat carry log-normalised values
in X and the counts in `.raw` (sometimes with the gene names lost on the way). `build_feature_map` turns an
eggNOG-mapper annotation of a proteome into the gene -> orthogroup table that `--feature-map` takes.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import numpy as np

_CHECK_CELLS = 500
DEFAULT_TAX_LEVEL = "33090"  # Viridiplantae, the level the dictionaries are built on


def _rows(X, n: int):
    sub = X[:n]
    return sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)


def _is_counts(X) -> bool:
    vals = X.data[: 200_000] if hasattr(X, "data") else _rows(X, _CHECK_CELLS).ravel()
    vals = np.asarray(vals, dtype=float)
    return bool(vals.size == 0 or (vals.min() >= 0 and np.allclose(vals, np.round(vals))))


def ensure_counts(adata: Any) -> tuple[Any, str]:
    """Return an AnnData whose X holds raw counts, and a note on what was done.

    Order: X already counts -> unchanged; a 'counts' layer -> used by the scoring as before; `.raw` with counts ->
    copied into X when its genes are X's genes (by name, or, when the names were lost in a conversion, by the same
    shape and the same non-zero pattern as X). Otherwise X is used as it is and the note says so."""
    import anndata as ad

    if _is_counts(adata.X):
        return adata, "X holds raw counts"
    if "counts" in adata.layers:
        return adata, "X is not raw counts; the 'counts' layer is used"
    raw = adata.raw
    if raw is not None and raw.shape == adata.shape and _is_counts(raw.X):
        same_names = list(map(str, raw.var_names)) == list(map(str, adata.var_names))
        n = min(_CHECK_CELLS, adata.n_obs)
        same_pattern = np.array_equal(_rows(raw.X, n) > 0, _rows(adata.X, n) > 0)
        if same_names or same_pattern:
            out = ad.AnnData(X=raw.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
            how = "gene names" if same_names else "matching non-zero pattern (gene names in .raw were missing)"
            return out, f"X is not raw counts; the counts in .raw are used (matched to X by {how})"
    return adata, ("WARNING X does not hold raw counts and no counts were found (.raw / 'counts' layer); scores are "
                   "computed on X as given, which can shift the calls - supply raw counts where possible")


def _select_og(field: str, levels: list[str]) -> str | None:
    by_tax: dict[str, str] = {}
    for ent in (field or "").split(","):
        if "@" not in ent:
            continue
        ogid, rest = ent.strip().split("@", 1)
        taxid = rest.split("|", 1)[0].strip()
        if ogid and taxid.isdigit():
            by_tax.setdefault(taxid, ogid)
    return next((by_tax[lvl] for lvl in levels if lvl in by_tax), None)


def build_feature_map(annotations: str | Path, species: str, out: str | Path,
                      strip_regex: str | None = r"\.\d+$", tax_level: str = DEFAULT_TAX_LEVEL) -> dict:
    """eggNOG-mapper `.emapper.annotations` -> gene -> orthogroup TSV (one orthogroup per gene, the best-scoring
    protein wins when isoforms disagree). `strip_regex` turns protein / transcript ids into gene ids."""
    strip = re.compile(strip_regex) if strip_regex else None
    best: dict[str, tuple] = {}
    n_rows = n_no_og = 0
    header = None
    with open(annotations, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            n_rows += 1
            query = (row.get("query") or "").strip()
            og = _select_og(row.get("eggNOG_OGs", ""), [tax_level])
            if not query or og is None:
                n_no_og += 1
                continue
            try:
                score = float(row.get("score", ""))
                conf = f"{score:.6g}"
            except ValueError:
                score, conf = float("-inf"), ""
            gene = strip.sub("", query) if strip else query
            feature = f"PCIP_OG:eggnog:{og}"
            cand = (score, feature, query, conf)
            prev = best.get(gene)
            if prev is None or cand[:3] > prev[:3]:
                best[gene] = cand
    if header is None:
        raise SystemExit(f"build-map: {annotations} does not look like an eggNOG-mapper .annotations file")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "matrix_gene_id", "canonical_id", "feature_id", "weight", "confidence", "source"])
        for gene, (_s, feature, query, conf) in best.items():
            w.writerow([species, gene, query, feature, "1.0", conf, "eggnog_viridiplantae"])
    return {"annotation_rows": n_rows, "rows_without_orthogroup": n_no_og, "genes_mapped": len(best),
            "orthogroups": len({v[1] for v in best.values()}), "out": str(out)}
