"""Gene -> cross-species orthogroup mapping and aggregation.

Maps a species gene-level matrix onto the cross-species orthogroup layer the dictionaries are written in. The
groups come from the ``gene_to_feature.tsv`` table (one row per gene -> orthogroup membership), columns:
    species, matrix_gene_id, canonical_id, feature_id, weight, confidence, source

  * Genes are joined by id (``var['gene_id']`` or var_names against ``matrix_gene_id``), never by position.
  * Per-gene ``var['orthomap_status']`` records mapped / ambiguous / unmapped; the alignment audit reports the
    mapping rate.
  * The input AnnData is never modified; aggregation returns a new AnnData.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

AGGREGATION_METHODS = ("sum", "max", "mean", "weighted")

FEATURE_MAP_COLUMNS = (
    "species",
    "matrix_gene_id",
    "canonical_id",
    "feature_id",
    "weight",
    "confidence",
    "source",
)

# orthogroup-layer var fields
OG_LAYER_VAR_FIELDS = (
    "orthogroup_id",
    "member_gene_ids",
    "member_species",
    "aggregation_method",
    "mapping_confidence",
)


# --- feature-map contract ----------------------------------------------------


@dataclass
class FeatureMap:
    """A loaded, validated gene->feature contract for one species (or all).

    ``by_gene[matrix_gene_id] = [(feature_id, weight, confidence), ...]``
    (a gene may belong to >1 feature for soft / many-to-many assignments).
    """

    by_gene: dict[str, list[tuple[str, float, float]]]
    source: str
    species: str | None = None

    def feature_ids(self) -> list[str]:
        out: set[str] = set()
        for members in self.by_gene.values():
            for fid, _w, _c in members:
                out.add(fid)
        return sorted(out)


def load_feature_map(path: str | Path, species: str | None = None) -> FeatureMap:
    """Load and validate the gene_to_feature.tsv contract.

    If ``species`` is given, only rows for that species are kept (the join is
    always within-species). ``weight`` defaults to 1.0, ``confidence`` to NaN
    when absent/blank.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"feature map not found: {path}")

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in ("species", "matrix_gene_id", "feature_id") if c not in header]
        if missing:
            raise ValueError(f"feature map {path}: missing required column(s) {missing}")

        by_gene: dict[str, list[tuple[str, float, float]]] = {}
        sources: set[str] = set()
        for row in reader:
            if species is not None and (row.get("species") or "").strip() != species:
                continue
            gid = (row.get("matrix_gene_id") or "").strip()
            fid = (row.get("feature_id") or "").strip()
            if not gid or not fid:
                continue
            weight = _to_float(row.get("weight"), default=1.0)
            conf = _to_float(row.get("confidence"), default=float("nan"))
            by_gene.setdefault(gid, []).append((fid, weight, conf))
            src = (row.get("source") or "").strip()
            if src:
                sources.add(src)

    return FeatureMap(by_gene=by_gene, source=",".join(sorted(sources)), species=species)


def _to_float(value: Any, default: float) -> float:
    if value is None:
        return default
    s = str(value).strip()
    if s == "" or s.lower() in ("na", "nan", "none"):
        return default
    try:
        return float(s)
    except ValueError:
        return default


# --- gene-id normalization (matrix_gene_id -> canonical_id) -------------------


@dataclass
class NormalizationReport:
    n_total: int = 0
    n_changed: int = 0
    n_ambiguous: int = 0
    ambiguous_ids: list[str] = field(default_factory=list)


def _load_alias_table(path: str | Path) -> dict[str, list[str]]:
    """Load a from_id<TAB>to_id crosswalk (many-to-one allowed)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"alias table not found: {path}")
    out: dict[str, list[str]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        # tolerate a header row named from_id/to_id, else treat as data
        if header and header[:2] != ["from_id", "to_id"]:
            if len(header) >= 2:
                out.setdefault(header[0], []).append(header[1])
        for row in reader:
            if len(row) >= 2 and row[0]:
                out.setdefault(row[0], []).append(row[1])
    return out


def normalize_gene_ids(
    gene_ids: list[str],
    species: str,
    rules: list[dict[str, Any]] | None,
    *,
    ambiguous_policy: str = "drop_with_flag",
) -> tuple[list[str | None], list[str], NormalizationReport]:
    """Apply per-species normalization rules to matrix gene ids.

    Returns (canonical_ids, statuses, report) where:
      canonical_ids[i] is the normalized id (or None if dropped as ambiguous),
      statuses[i] in {"identity","normalized","ambiguous"}.
    """
    rules = rules or []
    alias_cache: dict[str, dict[str, list[str]]] = {}
    canonical: list[str | None] = []
    statuses: list[str] = []
    report = NormalizationReport(n_total=len(gene_ids))

    for gid in gene_ids:
        cur = gid
        changed = False
        ambiguous = False
        for rule in rules:
            rtype = rule.get("type")
            if rtype == "strip_regex":
                new = re.sub(rule["pattern"], "", cur)
                if new != cur:
                    cur, changed = new, True
            elif rtype == "alias_table":
                table = alias_cache.get(rule["path"])
                if table is None:
                    table = _load_alias_table(rule["path"])
                    alias_cache[rule["path"]] = table
                hits = table.get(cur)
                if hits:
                    uniq = list(dict.fromkeys(hits))  # order-preserving dedupe
                    if len(uniq) > 1:
                        if ambiguous_policy == "keep_first":
                            # deterministically take the first target (not ambiguous)
                            if uniq[0] != cur:
                                cur, changed = uniq[0], True
                        else:
                            ambiguous = True
                            break
                    elif uniq[0] != cur:
                        cur, changed = uniq[0], True
        if ambiguous:
            report.n_ambiguous += 1
            if len(report.ambiguous_ids) < 100:
                report.ambiguous_ids.append(gid)
            if ambiguous_policy == "error":
                raise ValueError(f"{species}: ambiguous gene id '{gid}' under alias table")
            # drop_with_flag (default): record as ambiguous, drop the canonical id
            canonical.append(None)
            statuses.append("ambiguous")
            continue
        canonical.append(cur)
        if changed:
            report.n_changed += 1
            statuses.append("normalized")
        else:
            statuses.append("identity")
    return canonical, statuses, report


# --- annotate + strict alignment audit ---------------------------------------


@dataclass
class AlignmentAudit:
    dataset_id: str
    species: str
    n_genes: int
    n_mapped: int
    n_ambiguous: int
    n_unmapped: int
    match_rate: float
    expressed_match_rate: float
    n_features: int
    passed: bool
    layer_used: str
    expressed_gate_applicable: bool = True
    alignment_status: str = "ok"   # "ok" | "ok_degraded" (passed on gene-match, no usable counts)
    top_unmapped_expressed: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "dataset_id": self.dataset_id,
            "species": self.species,
            "n_genes": self.n_genes,
            "n_mapped": self.n_mapped,
            "n_ambiguous": self.n_ambiguous,
            "n_unmapped": self.n_unmapped,
            "match_rate": self.match_rate,
            "expressed_match_rate": self.expressed_match_rate,
            "expressed_gate_applicable": self.expressed_gate_applicable,
            "alignment_status": self.alignment_status,
            "n_features": self.n_features,
            "passed": self.passed,
            "layer_used": self.layer_used,
            "top_unmapped_expressed": self.top_unmapped_expressed,
        }


def annotate_and_audit(
    adata: "Any",
    feature_map: FeatureMap,
    *,
    dataset_id: str,
    species: str,
    normalization_rules: list[dict[str, Any]] | None = None,
    ambiguous_policy: str = "drop_with_flag",
    layer_priority: tuple[str, ...] = ("counts", "normalized"),
    min_match_rate: float = 0.30,
    min_expressed_match_rate: float = 0.50,
    expression_gate_policy: str = "require",
    raw_counts_available: bool | None = None,
    report_top_unmapped: int = 50,
) -> tuple["Any", AlignmentAudit]:
    """Annotate var with canonical_id/orthogroup_id/orthomap_status and run the
    strict alignment audit. Returns (annotated_copy, audit).

    The returned object is a copy; the input adata is never modified.
    """
    import numpy as np

    out = adata.copy()
    gene_ids = out.var["gene_id"].astype(str).tolist() if "gene_id" in out.var.columns \
        else list(out.var_names)

    canonical, norm_status, _nrep = normalize_gene_ids(
        gene_ids, species, normalization_rules, ambiguous_policy=ambiguous_policy
    )

    primary_feature: list[str] = []
    primary_conf: list[float] = []
    status: list[str] = []
    for raw_id, cid, nstat in zip(gene_ids, canonical, norm_status):
        if nstat == "ambiguous" and cid is None:
            primary_feature.append("")
            primary_conf.append(float("nan"))
            status.append("ambiguous")
            continue
        # the map is keyed by matrix gene id; fall back to the normalized id
        members = feature_map.by_gene.get(raw_id, [])
        if not members and cid and cid != raw_id:
            members = feature_map.by_gene.get(cid, [])
        if not members:
            primary_feature.append("")
            primary_conf.append(float("nan"))
            status.append("unmapped")
            continue
        # primary feature = highest weight, then highest confidence, then lexical
        best = sorted(
            members,
            key=lambda m: (m[1], (m[2] if m[2] == m[2] else -1.0), m[0]),
            reverse=True,
        )[0]
        primary_feature.append(best[0])
        primary_conf.append(best[2])
        status.append("mapped")

    out.var["canonical_id"] = [c if c is not None else "" for c in canonical]
    out.var["orthogroup_id"] = primary_feature
    out.var["orthogroup_confidence"] = primary_conf
    out.var["orthomap_status"] = status

    # expressed-weighted match: how much of total expressed signal is on mapped genes
    mapped_mask = np.array([s == "mapped" for s in status])
    n_mapped = int(mapped_mask.sum())
    n_amb = int(sum(1 for s in status if s == "ambiguous"))
    n_genes = len(status)
    n_unmapped = n_genes - n_mapped - n_amb
    match_rate = n_mapped / n_genes if n_genes else 0.0

    layer_used = next((l for l in layer_priority if l in out.layers), "X")
    X = out.layers[layer_used] if layer_used in out.layers else out.X
    gene_totals = np.asarray(X.sum(axis=0)).ravel()
    total = float(gene_totals.sum())
    expressed_match = float(gene_totals[mapped_mask].sum()) / total if total > 0 else 0.0

    # expressed_match needs counts: on centred data the column sums are ~0. With 'skip_if_no_counts' the
    # gate is skipped and the audit passes on gene match alone (status 'ok_degraded').
    if raw_counts_available is None:
        rc = out.uns.get("raw_counts_available")
        raw_counts_available = (str(rc).strip().lower() == "true") if rc is not None else None
    # counts declared present but summing to zero are empty data and still fail the gate
    no_usable_counts = (raw_counts_available is False) or (raw_counts_available is None and total <= 0)
    expressed_gate_applicable = not (expression_gate_policy == "skip_if_no_counts"
                                     and no_usable_counts)

    # most highly expressed unmapped genes
    unmapped_mask = np.array([s == "unmapped" for s in status])
    top_unmapped: list[str] = []
    if unmapped_mask.any():
        idx = np.where(unmapped_mask)[0]
        order = idx[np.argsort(-gene_totals[idx])][:report_top_unmapped]
        top_unmapped = [gene_ids[i] for i in order]

    n_features = len({f for f in primary_feature if f})
    if expressed_gate_applicable:
        passed = match_rate >= min_match_rate and expressed_match >= min_expressed_match_rate
        alignment_status = "ok" if passed else "failed_alignment"
    else:
        passed = match_rate >= min_match_rate   # gene-match only (no usable counts)
        alignment_status = "ok_degraded" if passed else "failed_alignment"

    audit = AlignmentAudit(
        dataset_id=dataset_id, species=species, n_genes=n_genes, n_mapped=n_mapped,
        n_ambiguous=n_amb, n_unmapped=n_unmapped, match_rate=match_rate,
        expressed_match_rate=expressed_match, n_features=n_features, passed=passed,
        layer_used=layer_used, expressed_gate_applicable=expressed_gate_applicable,
        alignment_status=alignment_status, top_unmapped_expressed=top_unmapped,
    )
    out.uns["orthomap_audit"] = audit.to_dict()
    return out, audit


# --- aggregation to the feature layer ----------------------------------------


def aggregate_to_features(
    adata: "Any",
    feature_map: FeatureMap,
    *,
    method: str = "sum",
    species: str = "",
    layer_priority: tuple[str, ...] = ("counts", "normalized"),
) -> "Any":
    """Aggregate a gene-layer AnnData onto the feature layer (cells x features).

    Uses the var['canonical_id'] produced by ``annotate_and_audit`` (falls back
    to gene_id). sum/mean/weighted are sparse incidence-matrix products; max is a
    grouped reduction. Returns a new AnnData with the orthogroup-layer var schema.
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    if method not in AGGREGATION_METHODS:
        raise ValueError(f"unknown aggregation method '{method}' (allowed: {AGGREGATION_METHODS})")

    layer_used = next((l for l in layer_priority if l in adata.layers), "X")
    X = adata.layers[layer_used] if layer_used in adata.layers else adata.X
    X = sp.csc_matrix(X)  # cells x genes; column slicing for max

    canon = adata.var["canonical_id"].astype(str).tolist() if "canonical_id" in adata.var.columns \
        else adata.var["gene_id"].astype(str).tolist()

    # Build gene-index -> list of (feature, weight) using the map.
    features = feature_map.feature_ids()
    fidx = {f: j for j, f in enumerate(features)}
    n_cells, n_genes = X.shape
    n_feat = len(features)

    # member gene id / species bookkeeping for var schema
    members_of: dict[str, list[str]] = {f: [] for f in features}

    if method == "max":
        result = sp.lil_matrix((n_cells, n_feat), dtype=np.float32)
        Xcsc = X.tocsc()
        # accumulate per feature the elementwise max across its member gene columns
        feat_cols: dict[int, list[int]] = {}
        for gi, cid in enumerate(canon):
            for f, _w, _c in feature_map.by_gene.get(cid, []):
                feat_cols.setdefault(fidx[f], []).append(gi)
                members_of[f].append(adata.var_names[gi])
        for fj, cols in feat_cols.items():
            sub = Xcsc[:, cols]
            mx = sub.max(axis=1)  # sparse max over columns -> column vector
            result[:, fj] = mx
        result = result.tocsr()
    else:
        # sum / mean / weighted via a sparse genes x features incidence matrix
        rows, cols, vals = [], [], []
        feat_gene_count = np.zeros(n_feat, dtype=np.float64)
        for gi, cid in enumerate(canon):
            for f, w, _c in feature_map.by_gene.get(cid, []):
                fj = fidx[f]
                rows.append(gi)
                cols.append(fj)
                vals.append(w if method == "weighted" else 1.0)
                feat_gene_count[fj] += 1.0
                members_of[f].append(adata.var_names[gi])
        incidence = sp.csr_matrix((vals, (rows, cols)), shape=(n_genes, n_feat), dtype=np.float32)
        result = (X.tocsr() @ incidence)  # cells x features
        if method == "mean":
            with np.errstate(divide="ignore", invalid="ignore"):
                inv = np.where(feat_gene_count > 0, 1.0 / feat_gene_count, 0.0)
            result = result @ sp.diags(inv.astype(np.float32))
        result = sp.csr_matrix(result)

    out = ad.AnnData(X=result.astype(np.float32))
    out.obs = adata.obs.copy()
    out.obs_names = adata.obs_names
    out.var_names = list(features)
    # orthogroup-layer var fields
    member_species = species if species else (
        str(adata.obs["species"].iloc[0]) if "species" in adata.obs.columns else ""
    )
    out.var["orthogroup_id"] = list(features)
    out.var["member_gene_ids"] = [",".join(members_of[f]) for f in features]
    out.var["member_species"] = member_species
    out.var["aggregation_method"] = method
    out.var["mapping_confidence"] = np.nan
    out.layers[f"orthogroup_{method}"] = out.X.copy()

    out.uns["orthomap"] = {
        "method": method,
        "layer_used": layer_used,
        "n_features": n_feat,
        "feature_map_source": feature_map.source,
    }
    if "orthomap_audit" in adata.uns:
        out.uns["orthomap_audit"] = adata.uns["orthomap_audit"]
    return out

