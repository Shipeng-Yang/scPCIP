"""scPCIP command-line interface.

    scpcip predict --input query.h5ad --species "Arabidopsis thaliana" --out results/ [--organ-key organ] [--no-report]

One predicted cell type per cell (with tissue system, confidence columns and program scores) from the shipped
root / leaf / flower / seed dictionaries.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from . import resources

ORGAN_CONFIGS = {"root": "predict_root.yaml", "leaf": "predict_leaf.yaml",
                 "flower": "predict_flower.yaml", "seed": "predict_seed.yaml"}
DEFAULT_ORGAN = "root"
MIN_ORGAN_PART_CELLS = 20  # each organ part is predicted on its own; the refinement needs enough cells
UNSUPPORTED_ORGANS = {"stem": "the stem dictionary is below the single-marker baseline on held-out species and is "
                              "not shipped"}


def _package_versions() -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    versions: dict[str, str] = {"scpcip": __version__}
    for dist in ("numpy", "scipy", "pandas", "anndata", "pyyaml", "scikit-learn"):
        try:
            versions[dist] = version(dist)
        except PackageNotFoundError:
            versions[dist] = "not_installed"
    return versions


def _dump_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False)


def _dataset_stem(path: Path) -> str:
    name = path.name
    for suffix in (".sum.h5ad", ".h5ad"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _sha256_file(path: str | Path | None) -> str | None:
    if not path:
        return None
    import hashlib

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_h5ad(path: Path):
    import anndata as ad

    return ad.read_h5ad(str(path))


def _read_organ_map(path: str) -> dict[str, str]:
    """Two tab-separated columns: sample value -> organ ('#' lines and blank lines ignored)."""
    out = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            raise SystemExit(f"predict: --organ-map line {n}: expected 'sample<TAB>organ', got {line!r}")
        out[parts[0].strip()] = parts[1].strip()
    return out


def _resolve_config(args: argparse.Namespace) -> tuple[Path, dict]:
    """--config (a file) or --organ (a shipped organ); paths inside a shipped config are package resources."""
    import yaml

    if args.config is not None:
        path = Path(args.config)
        if not path.exists():
            raise SystemExit(f"predict: config not found: {path}")
    else:
        organ = (args.organ or DEFAULT_ORGAN).lower()
        if organ in UNSUPPORTED_ORGANS:
            raise SystemExit(f"predict: organ {organ!r} is not supported: {UNSUPPORTED_ORGANS[organ]}")
        if organ not in ORGAN_CONFIGS:
            raise SystemExit(f"predict: unknown organ {organ!r}; supported: {', '.join(ORGAN_CONFIGS)}")
        path = resources.config_path(ORGAN_CONFIGS[organ])
    cfg = yaml.safe_load(path.read_text()) or {}
    cfg["dictionary"] = str(resources.resolve(cfg["dictionary"], base=path.parent))
    return path, cfg


def _predict_by_organ(args: argparse.Namespace) -> int:
    """Split the query by the user-given organ of each cell, predict each part with that organ's shipped config and
    merge (input cell order kept, `organ` column added). Fails on a missing column, an unmapped value, an unknown
    organ, an unsupported organ (stem) or an organ part with too few cells."""
    import copy
    import shutil
    import tempfile
    import pandas as pd

    if args.config is not None or args.organ is not None:
        raise SystemExit("predict: --config / --organ cannot be combined with --organ-key")
    in_path = Path(args.input)
    adata = _read_h5ad(in_path)
    if args.organ_key not in adata.obs.columns:
        raise SystemExit(f"predict: --organ-key column {args.organ_key!r} not in the input obs "
                         f"(columns: {', '.join(map(str, adata.obs.columns))})")
    if not adata.obs_names.is_unique:
        raise SystemExit("predict: --organ-key needs unique cell ids (obs_names)")
    values = adata.obs[args.organ_key].astype(str).to_numpy()
    if args.organ_map:
        omap = _read_organ_map(args.organ_map)
        missing = sorted(set(values) - set(omap))
        if missing:
            raise SystemExit(f"predict: --organ-map has no organ for {args.organ_key} value(s): {', '.join(missing)}")
        organs = [omap[v] for v in values]
    else:
        organs = list(values)
    organs = pd.Series([str(o).strip().lower() for o in organs], index=adata.obs_names)
    bad = sorted(set(organs) - set(ORGAN_CONFIGS))
    for o in bad:
        if o in UNSUPPORTED_ORGANS:
            raise SystemExit(f"predict: organ {o!r} is not supported: {UNSUPPORTED_ORGANS[o]}")
    if bad:
        raise SystemExit(f"predict: unknown organ(s) {', '.join(bad)}; supported: {', '.join(ORGAN_CONFIGS)}")
    small = {o: int(n) for o, n in organs.value_counts().items() if n < MIN_ORGAN_PART_CELLS}
    if small:
        raise SystemExit("predict: too few cells to predict an organ part (need >= "
                         f"{MIN_ORGAN_PART_CELLS}): " + ", ".join(f"{o} {n} cells" for o, n in sorted(small.items())))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    dsid = args.dataset_id or _dataset_stem(in_path)
    runs, preds, scores = {}, [], []
    tmp = Path(tempfile.mkdtemp(prefix=".organ_split_", dir=out_dir))
    try:
        for organ in [o for o in ORGAN_CONFIGS if (organs == o).any()]:
            part = adata[(organs == organ).to_numpy()].copy()
            part_path = tmp / f"{dsid}__{organ}.h5ad"
            part.write_h5ad(part_path)
            sub = copy.copy(args)
            sub.input, sub.config, sub.organ, sub.organ_key, sub.organ_map = str(part_path), None, organ, None, None
            sub.dataset_id = f"{dsid}__{organ}"
            if _predict_run(sub) != 0:
                raise SystemExit(f"predict: the {organ} part failed")
            pr = pd.read_csv(out_dir / f"{sub.dataset_id}.predictions.tsv", sep="\t", keep_default_na=False)
            pr.insert(2, "organ", organ)
            preds.append(pr)
            sc = pd.read_csv(out_dir / f"{sub.dataset_id}.type_scores.tsv", sep="\t", index_col=0)
            sc.insert(0, "organ", organ)
            scores.append(sc)
            runs[organ] = {"n_cells": int(part.n_obs), "config": ORGAN_CONFIGS[organ], "dataset_id": sub.dataset_id,
                           "report": str(out_dir / f"{sub.dataset_id}.predict_report.json")}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    order = list(map(str, adata.obs_names))
    merged = pd.concat(preds, ignore_index=True).set_index("cell_id").loc[order].reset_index()
    merged["dataset_id"] = dsid
    merged.to_csv(out_dir / f"{dsid}.predictions.tsv", sep="\t", index=False)
    pd.concat(scores).loc[order].to_csv(out_dir / f"{dsid}.type_scores.tsv", sep="\t", index_label="cell_id")
    report = {"command": "predict", "mode": "by_organ", "argv": sys.argv, "input": str(in_path),
              "input_sha256": _sha256_file(in_path), "dataset_id": dsid, "species": args.species,
              "organ_key": args.organ_key, "organ_map": args.organ_map,
              "organ_map_sha256": _sha256_file(args.organ_map), "runs": runs, "n_cells": len(order),
              "note": "each organ part was scored on its own (score correction and confidence quartiles within it); "
                      "per-organ outputs are <dataset_id>__<organ>.*",
              "finished_utc": datetime.now(timezone.utc).isoformat()}
    _dump_json(out_dir / f"{dsid}.predict_report.json", report)
    print(f"predict by organ done: {dsid} ({len(order)} cells; "
          + ", ".join(f"{o}={r['n_cells']}" for o, r in runs.items()) + ")")
    print(f"  predictions: {out_dir / f'{dsid}.predictions.tsv'}")
    return 0


def _cmd_predict(args: argparse.Namespace) -> int:
    # One BLAS / OpenMP thread: the self-training refinement otherwise changes a small fraction of calls with the
    # machine's thread count. The variables cover libraries not yet loaded; threadpool_limits covers loaded ones.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = "1"
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    from threadpoolctl import threadpool_limits

    if not args.og_level and not args.species:
        raise SystemExit("predict: --species is required unless --og-level is given")
    with threadpool_limits(limits=1):
        if args.organ_key:
            return _predict_by_organ(args)
        return _predict_run(args)


def _predict_run(args: argparse.Namespace) -> int:
    """gene -> orthogroup (unless --og-level) -> program scoring -> one predicted cell type per cell (argmax; seeded
    tie-break, flagged) with confidence columns. Never abstains."""
    import pandas as pd
    from . import predict as pred_mod
    from .program import load_dictionary

    config_path, cfg = _resolve_config(args)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    in_path = Path(args.input)
    dsid = args.dataset_id or _dataset_stem(in_path)
    dictionary = load_dictionary(cfg["dictionary"])
    seed = int(args.seed if args.seed is not None else cfg.get("seed", pred_mod.DEFAULT_SEED))
    tissue_map = {str(k): str(v) for k, v in (cfg.get("tissue_system") or {}).items()}
    if tissue_map:
        unmapped = sorted({str(p.get("cell_type_coarse")) for p in dictionary["programs"]} - set(tissue_map))
        if unmapped:
            raise SystemExit(f"predict: tissue_system in {config_path} has no entry for {', '.join(unmapped)}")

    adata = _read_h5ad(in_path)
    from .inputs import ensure_counts

    adata, counts_note = ensure_counts(adata)
    if counts_note.startswith("WARNING") or "not raw counts" in counts_note:
        print(f"predict: {counts_note}")
    report = {"command": "predict", "argv": sys.argv, "input": str(in_path),
              "input_sha256": _sha256_file(in_path), "dataset_id": dsid, "species": args.species,
              "config": str(config_path), "dictionary": Path(cfg["dictionary"]).name,
              "dictionary_sha256": _sha256_file(cfg["dictionary"]), "seed": seed,
              "og_level_input": bool(args.og_level), "input_counts": counts_note,
              "package_versions": _package_versions()}
    if args.og_level:
        og_adata = adata
    else:
        from . import orthology as ortho

        fmap_path = Path(args.feature_map) if args.feature_map else resources.feature_map_path()
        if not fmap_path.exists():
            raise SystemExit(f"predict: feature map not found: {fmap_path}")
        report["feature_map"] = fmap_path.name if not args.feature_map else str(fmap_path)
        report["feature_map_sha256"] = _sha256_file(fmap_path)
        if not args.feature_map and not os.environ.get("SCPCIP_FEATURE_MAP"):
            report["feature_map_revision"] = resources.FEATURE_MAP_REVISION
        fmap = ortho.load_feature_map(str(fmap_path), species=args.species)
        if not fmap.by_gene:
            raise SystemExit(f"predict: no gene->orthogroup mapping for species {args.species!r} in {fmap_path} - "
                             "check the species name (`scpcip species` lists the shipped ones) or supply a map for "
                             "it with --feature-map")
        annotated, audit = ortho.annotate_and_audit(adata, fmap, dataset_id=dsid, species=args.species,
                                                    min_match_rate=0.0, min_expressed_match_rate=0.0)
        if audit.match_rate <= 0:
            raise SystemExit(f"predict: none of the {adata.n_vars} input genes maps to an orthogroup for "
                             f"{args.species!r} - check the gene id namespace (e.g. {list(adata.var_names[:3])})")
        og_adata = ortho.aggregate_to_features(annotated, fmap, method=cfg.get("aggregate", "mean"),
                                               species=args.species)
        report["alignment_audit"] = audit.to_dict()
        report["og_mapping_rate"] = float(audit.match_rate)
        warn_below = float(cfg.get("warn_mapping_rate_below", 0.5))
        if audit.match_rate < warn_below:
            print(f"predict: WARNING gene->orthogroup mapping rate {audit.match_rate:.3f} < {warn_below} - "
                  f"check --species and the gene id namespace; predictions may be unreliable")
    preds, scores, res = pred_mod.predict_og_adata(og_adata, dictionary, method=cfg.get("method", "weighted_rank"),
                                                   seed=seed)
    report_preds = preds
    hierarchical = (bool(dictionary.get("fine_programs")) or bool(dictionary.get("state_programs"))
                    or dictionary.get("hierarchical") is True)
    report["hierarchical"] = hierarchical
    if hierarchical:
        # per-query score correction + label-free refinement at the coarse level; fine types are used internally by
        # the hierarchy but not reported as the cell type
        from . import hierarchy as hi
        from .report import margin_quartile

        corr = str(cfg.get("score_correction", dictionary.get("default_score_correction", "bg_center")))
        refine = str(cfg.get("refine", "none"))
        hp, scores = hi.predict_hierarchy(og_adata, dictionary, score_correction=corr,
                                          refine=refine, vote_resolution=float(cfg.get("vote_resolution", 2.0)),
                                          vote_min_purity=float(cfg.get("vote_min_purity", 0.5)),
                                          refine_growth_cap=(None if cfg.get("refine_growth_cap") is None
                                                             else float(cfg["refine_growth_cap"])),
                                          return_scores=True, method=cfg.get("method", "weighted_rank"), seed=seed,
                                          return_fine_scores=bool(dictionary.get("reference_subtypes")))
        # exported type_scores = the corrected coarse scores the final labels and margins refer to
        preds = pd.DataFrame({"cell_id": hp["cell_id"].to_numpy(),
                              "predicted_cell_type": hp["predicted_coarse"].to_numpy(),
                              "predicted_coarse": hp["predicted_coarse"].to_numpy(),
                              "margin": hp["coarse_margin"].to_numpy(),
                              "tie_at_top": hp["coarse_tie_at_top"].to_numpy(),
                              "invalid_score": preds["invalid_score"].to_numpy()})
        report_preds = preds
        for c in [c for c in hp.columns if c.startswith("state_") and c.endswith("_score")]:
            preds[c.replace("_qc_initial", "")] = hp[c].to_numpy()  # developmental-state scores (never a cell type)
        # per-cell confidence = margin quartile within this query, with the accuracy measured for that quartile on
        # held-out species
        preds["confidence_quartile"] = margin_quartile(preds["margin"])
        qacc = (dictionary.get("reliability") or {}).get("margin_quartile_accuracy")
        preds["expected_accuracy"] = [qacc[q - 1] for q in preds["confidence_quartile"]] if qacc else float("nan")
        # optional secondary subtype inside listed coarse types (leaf: xylem / phloem inside vascular); the main call
        # is unchanged and every subtype score is reported
        ref_types = [str(t) for t in dictionary.get("reference_subtypes") or []]
        if ref_types:
            import numpy as np

            inside = hp["predicted_coarse"].isin(ref_types).to_numpy()
            preds["reference_subtype"] = np.where(inside, hp["predicted_fine"].astype(str).to_numpy(), "")
            preds["reference_subtype_margin"] = np.where(inside, hp["fine_margin"].to_numpy(dtype=float), np.nan)
            for t in ref_types:
                for col in [c for c in hp.columns if c.startswith(f"fine_score_{t}__")]:
                    preds[f"reference_subtype_score_{col[len(f'fine_score_{t}__'):]}"] = hp[col].to_numpy(dtype=float)
            report["reference_subtypes"] = ref_types
            report["reference_subtype_note"] = ((dictionary.get("reliability") or {}).get("reference_subtype")
                                                or "secondary, indicative only; not part of the main cell-type call")
        if "predicted_coarse_unrefined" in hp:
            preds["predicted_coarse_unrefined"] = hp["predicted_coarse_unrefined"].to_numpy()
            preds["margin_unrefined"] = hp["coarse_margin_unrefined"].to_numpy()
            preds["refine_changed"] = (hp["predicted_coarse"] != hp["predicted_coarse_unrefined"]).to_numpy()
        report.update({"score_correction": corr, "refine": refine,
                       "state_programs": [p["cell_type_coarse"] for p in dictionary.get("state_programs") or []],
                       "type_scores_are": f"coarse program scores after score_correction={corr}",
                       "confidence": "margin quartile within this query (1 = least confident); expected_accuracy = "
                                     "leave-one-species-out accuracy of that quartile"})
    preds.insert(1, "dataset_id", dsid)
    if tissue_map:
        preds.insert(list(preds.columns).index("predicted_cell_type") + 1, "tissue_system",
                     preds["predicted_cell_type"].astype(str).map(tissue_map).fillna("").to_numpy())
        report["tissue_system_map"] = tissue_map
        if cfg.get("tissue_system_reliability"):
            report["tissue_system_reliability"] = str(cfg["tissue_system_reliability"])
        report["tissue_system_counts"] = {str(k): int(v) for k, v in preds["tissue_system"].value_counts().items()}
    pred_path = out_dir / f"{dsid}.predictions.tsv"
    score_path = out_dir / f"{dsid}.type_scores.tsv"
    preds.to_csv(pred_path, sep="\t", index=False)
    scores.to_csv(score_path, sep="\t", index_label="cell_id")
    report.update({
        "n_cells": int(len(preds)),
        "predicted_counts": {str(k): int(v) for k, v in preds["predicted_cell_type"].value_counts().items()},
        "tie_rate": float(preds["tie_at_top"].mean()) if "tie_at_top" in preds else None,
        "invalid_score_cells": int(preds["invalid_score"].sum()),
        "program_coverage": {pid: float(c) for pid, c in zip(res.program_ids, res.coverage)},
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    })
    if args.report:
        from . import report as report_mod

        report["report_files"] = report_mod.build_report(og_adata, report_preds, scores, dictionary, out_dir, dsid, seed,
                                                         meta={k: report[k] for k in ("dataset_id", "species",
                                                               "dictionary", "seed", "input") if k in report},
                                                         full_preds=preds if hierarchical else None,
                                                         sample_key=args.sample_key, cluster_key=args.cluster_key)
        report["sample_consistency"] = report["report_files"].pop("sample_consistency_summary")
        if args.r_figures is not False:
            from . import rplots

            organ = None if args.config is not None else (args.organ or DEFAULT_ORGAN).lower()
            report["r_figures"] = rplots.make(adata, out_dir, dsid, organ, args.species)
    _dump_json(out_dir / f"{dsid}.predict_report.json", report)
    print(f"predict done: {dsid} ({len(preds)} cells, {len(dictionary['programs'])} programs)")
    print("  " + ", ".join(f"{k}={v}" for k, v in report["predicted_counts"].items()))
    if "og_mapping_rate" in report:
        print(f"  gene->orthogroup mapping rate: {report['og_mapping_rate']:.3f}")
    print(f"  predictions: {pred_path}")
    if args.report:
        print(f"  report:      {report['report_files']['html']}")
        cons = report["sample_consistency"]
        if cons["status"] in ("consistent", "check", "not_reproducible"):
            print(f"  reproducibility across samples ({cons['sample_key']}, clusters {cons['cluster_key']}): "
                  f"{cons['status']} (mean difference {cons['weighted_mean_tv']:.2f}, "
                  f"{cons['fraction_cells_assessed']:.0%} of cells assessed)")
            if cons["status"] == "not_reproducible":
                print("predict: WARNING the same clusters received different cell types in different samples - "
                      "do not use these calls without independent confirmation (see the report)")
        elif cons["status"] == "not_assessable":
            print(f"  reproducibility across samples: not assessable - {cons['reason']}")
        _print_r_figures(report.get("r_figures"), required=bool(args.r_figures))
    return 0


def _print_r_figures(res: dict | None, required: bool) -> None:
    if not res:
        return
    if res["status"] == "done":
        print(f"  R figures:   {res['dir']}")
    elif res["status"] == "failed" or required:
        print(f"predict: WARNING R figures not drawn - {res['reason']}")
    else:
        print(f"  R figures:   not drawn ({res['reason']})")


def _cmd_species(args: argparse.Namespace) -> int:
    """List the species covered by the shipped gene -> orthogroup map."""
    import csv

    path = Path(args.feature_map) if args.feature_map else resources.feature_map_path()
    counts: dict[str, int] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["species"]] = counts.get(row["species"], 0) + 1
    for sp in sorted(counts):
        print(f"{sp}\t{counts[sp]} genes")
    return 0


def _cmd_build_map(args: argparse.Namespace) -> int:
    """eggNOG-mapper annotation of a proteome -> gene -> orthogroup table for --feature-map."""
    from .inputs import build_feature_map

    info = build_feature_map(args.annotations, args.species, args.out,
                             strip_regex=None if args.no_strip else args.strip_regex)
    print(f"build-map: {info['genes_mapped']} genes -> {info['orthogroups']} orthogroups "
          f"({info['rows_without_orthogroup']} of {info['annotation_rows']} proteins had no Viridiplantae "
          f"orthogroup); written to {info['out']}")
    print(f'  use it with: scpcip predict --species "{args.species}" --feature-map {info["out"]} ...')
    return 0


def _cmd_plot(args: argparse.Namespace) -> int:
    """R figures (organ map, UMAPs, marker heatmap / violins) for a finished predict run."""
    from . import rplots
    from .inputs import ensure_counts

    exe, why = rplots.check()
    if exe is None:
        raise SystemExit(f"plot: {why}")
    in_path = Path(args.input)
    dsid = args.dataset_id or _dataset_stem(in_path)
    adata, _ = ensure_counts(_read_h5ad(in_path))
    res = rplots.make(adata, args.out, dsid, (args.organ or DEFAULT_ORGAN).lower(), args.species, exe=exe)
    if res["status"] != "done":
        raise SystemExit(f"plot: {res['reason']}")
    print(f"plot done: {res['dir']} ({len(res['files'])} files)")
    return 0


def _cmd_fetch(args: argparse.Namespace) -> int:
    """Download the gene -> orthogroup map into the local cache (done automatically on first use)."""
    print(resources.feature_map_path())
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scpcip", description="scPCIP: program-based plant single-cell annotation",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=resources.version_text(__version__))
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("predict", help="one cell type per cell (argmax over programs) + tissue system + confidence")
    p.add_argument("--input", required=True, help="query .h5ad (raw counts; gene ids in var_names)")
    p.add_argument("--species", default=None, help="query species, e.g. 'Arabidopsis thaliana' (see `scpcip species`)")
    p.add_argument("--organ", default=None, help="root (default) / leaf / flower / seed")
    p.add_argument("--organ-key", dest="organ_key", default=None,
                   help="obs column giving each cell's organ, or its sample when --organ-map is given: each organ "
                        "part is predicted with that organ's dictionary")
    p.add_argument("--organ-map", dest="organ_map", default=None,
                   help="two-column TSV: <value of --organ-key> <TAB> <organ>")
    p.add_argument("--config", default=None, help="a predict config YAML instead of a shipped organ")
    p.add_argument("--feature-map", dest="feature_map", default=None,
                   help="your own gene -> orthogroup TSV (columns species / matrix_gene_id / feature_id) for species "
                        "the shipped map does not cover")
    p.add_argument("--og-level", dest="og_level", action="store_true",
                   help="input var_names are already orthogroup ids (skip gene -> orthogroup mapping)")
    p.add_argument("--out", default="scpcip_out", help="output directory")
    p.add_argument("--dataset-id", dest="dataset_id", default=None, help="output file prefix (default: input name)")
    p.add_argument("--seed", type=int, default=None, help="tie-break seed (default: from the config)")
    p.add_argument("--report", action=argparse.BooleanOptionalAction, default=True,
                   help="write figures (PNG + PDF), marker / composition tables and a self-contained HTML report "
                        "(default; --no-report writes the prediction tables only)")
    p.add_argument("--sample-key", dest="sample_key", default=None,
                   help="obs column with each cell's sample / replicate for the reproducibility check in the report "
                        "(default: the first of sample, orig.ident, batch, library, replicate that is present)")
    p.add_argument("--cluster-key", dest="cluster_key", default=None,
                   help="obs column with clusters from an integration of the samples, for that check (default: "
                        "seurat_clusters / leiden / louvain / clusters / cluster if present; without one the check "
                        "is not run)")
    p.add_argument("--r-figures", dest="r_figures", action=argparse.BooleanOptionalAction, default=None,
                   help="also draw the R figures (organ map, UMAPs, marker heatmap) into <dataset_id>.scPCIP_figures/; "
                        "default: when Rscript with SCP, Seurat and ggPlantmap is available")
    p.set_defaults(func=_cmd_predict)

    p = sub.add_parser("plot", help="R figures (organ map, UMAPs, marker heatmap / violins) for a finished predict run")
    p.add_argument("--input", required=True, help="the query .h5ad given to predict")
    p.add_argument("--out", default="scpcip_out", help="the predict output directory")
    p.add_argument("--dataset-id", dest="dataset_id", default=None, help="as given to predict (default: input name)")
    p.add_argument("--organ", default=None, help="root (default) / leaf / flower / seed, as given to predict")
    p.add_argument("--species", default=None, help="query species (picks the grass or the dicot leaf drawing)")
    p.set_defaults(func=_cmd_plot)

    p = sub.add_parser("species", help="list the species covered by the shipped gene -> orthogroup map")
    p.add_argument("--feature-map", dest="feature_map", default=None)
    p.set_defaults(func=_cmd_species)

    p = sub.add_parser("build-map", help="gene -> orthogroup table for a species not in the shipped map, from an "
                                         "eggNOG-mapper annotation of its proteins")
    p.add_argument("--annotations", required=True, help="<name>.emapper.annotations from eggNOG-mapper")
    p.add_argument("--species", required=True, help="species name to use with --species in predict")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--strip-regex", dest="strip_regex", default=r"\.\d+$",
                   help="regex removed from protein ids to get the gene ids used in the count matrix "
                        "(default: a trailing isoform number such as .1)")
    p.add_argument("--no-strip", dest="no_strip", action="store_true", help="protein ids are already gene ids")
    p.set_defaults(func=_cmd_build_map)

    p = sub.add_parser("fetch", help="download the gene -> orthogroup map now (otherwise fetched on first use)")
    p.set_defaults(func=_cmd_fetch)
    return parser


def main(argv: list[str] | None = None) -> int:
    from . import banner
    if banner.enabled() and sys.stderr.isatty():
        banner.show()
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
