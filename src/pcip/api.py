"""Python interface: annotate an in-memory AnnData.

    import scanpy as sc, pcip
    adata = sc.read_h5ad("root.h5ad")
    pcip.annotate(adata, species="Oryza sativa", organ="root")
    adata.obs["scpcip_predicted_cell_type"]

The data are written to a temporary .h5ad and passed through exactly the same code path as `scpcip predict`, so the
results are identical to the command line.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
_DATASET_ID = "query"


def annotate(adata: Any, species: str | None = None, organ: str | None = None, *, organ_key: str | None = None,
             organ_map: str | None = None, feature_map: str | None = None, og_level: bool = False,
             seed: int | None = None, out: str | Path | None = None, report: bool | None = None,
             sample_key: str | None = None, cluster_key: str | None = None,
             key_added: str = "scpcip", copy: bool = False) -> Any:
    """Predict one cell type per cell and store the result in ``adata``.

    Parameters mirror ``scpcip predict``: ``organ`` is root (default) / leaf / flower / seed; ``organ_key`` names an
    obs column with each cell's organ (or sample, with ``organ_map``). With ``out`` the output files are kept and the
    figures and HTML report are written there (``report=False`` skips them); without it nothing is left on disk.

    Adds ``obs['<key_added>_<column>']`` for every column of ``predictions.tsv`` (predicted_cell_type, tissue_system,
    margin, confidence_quartile, expected_accuracy, ...), ``obsm['<key_added>_type_scores']`` (cells x types) and
    ``uns['<key_added>']`` (the run report). Returns the annotated copy when ``copy=True``, otherwise None.
    """
    import pandas as pd

    from .cli import main

    if report is None:
        report = out is not None
    if report and out is None:
        raise ValueError("annotate: report=True needs out= (a directory for the figures and HTML report)")
    if not adata.obs_names.is_unique:
        raise ValueError("annotate: cell ids (obs_names) must be unique")
    if copy:
        adata = adata.copy()

    saved_env = {v: os.environ.get(v) for v in _THREAD_VARS}
    with tempfile.TemporaryDirectory(prefix="scpcip_") as tmp:
        in_path = Path(tmp) / f"{_DATASET_ID}.h5ad"
        adata.write_h5ad(in_path)
        out_dir = Path(out) if out is not None else Path(tmp) / "out"
        argv = ["predict", "--input", str(in_path), "--out", str(out_dir), "--dataset-id", _DATASET_ID]
        for flag, value in (("--species", species), ("--organ", organ), ("--organ-key", organ_key),
                            ("--organ-map", organ_map), ("--feature-map", feature_map), ("--seed", seed),
                            ("--sample-key", sample_key), ("--cluster-key", cluster_key)):
            if value is not None:
                argv += [flag, str(value)]
        if og_level:
            argv.append("--og-level")
        argv.append("--report" if report else "--no-report")
        try:
            main(argv)
        except SystemExit as exc:
            if exc.code in (0, None):
                raise
            raise ValueError(f"annotate: {exc.code}") from None
        finally:
            for var, value in saved_env.items():
                if value is None:
                    os.environ.pop(var, None)
                else:
                    os.environ[var] = value
        preds = pd.read_csv(out_dir / f"{_DATASET_ID}.predictions.tsv", sep="\t")
        scores = pd.read_csv(out_dir / f"{_DATASET_ID}.type_scores.tsv", sep="\t", index_col=0)
        report_path = out_dir / f"{_DATASET_ID}.predict_report.json"
        run_report = json.loads(report_path.read_text())
        run_report["input"] = "AnnData passed to pcip.annotate"  # the temporary .h5ad is deleted on return
        if out is not None:
            report_path.write_text(json.dumps(run_report, indent=1))

    cells = adata.obs_names.astype(str)
    preds = preds.set_index(preds["cell_id"].astype(str)).loc[cells]
    for col in preds.columns.drop(["cell_id", "dataset_id"]):
        adata.obs[f"{key_added}_{col}"] = preds[col].to_numpy()
    scores.index = scores.index.astype(str)
    scores = scores.loc[cells]
    scores.index = adata.obs_names
    adata.obsm[f"{key_added}_type_scores"] = scores
    adata.uns[key_added] = run_report
    return adata if copy else None
