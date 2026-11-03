"""Optional R figures: organ map (ggPlantmap), UMAPs, marker heatmap and violins (SCP).

Drawn from the tables a predict run with the report writes ({id}.predictions.tsv, {id}.umap.tsv, {id}.markers.tsv)
plus the query expression of the marker member genes, into ``{out}/{id}.scPCIP_figures/``. Needs Rscript with the
R packages SCP, Seurat and ggPlantmap (all preinstalled in the scPCIP container image); set SCPCIP_RSCRIPT to choose
the Rscript. Nothing here changes predictions; when R is missing the run goes on without these figures.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from importlib import resources as _ir
from pathlib import Path
from typing import Any

import numpy as np

R_PACKAGES = ("SCP", "Seurat", "ggPlantmap")
TIMEOUT_S = 1800
# grass leaves are drawn on the rice leaf cross-section; every other leaf on the Arabidopsis one
GRASS_GENERA = {"Oryza", "Zea", "Sorghum", "Triticum", "Hordeum", "Setaria", "Brachypodium", "Panicum", "Saccharum",
                "Avena", "Secale", "Eleusine", "Pennisetum", "Cenchrus", "Lolium", "Festuca", "Miscanthus",
                "Phyllostachys", "Eragrostis", "Digitaria", "Leersia", "Aegilops"}


def script_path() -> Path:
    return Path(str(_ir.files("pcip") / "data" / "R" / "scpcip_figures.R"))


def rscript() -> str | None:
    return os.environ.get("SCPCIP_RSCRIPT") or shutil.which("Rscript")


def check() -> tuple[str | None, str]:
    """(Rscript, "") when the R figures can be drawn, else (None, reason)."""
    exe = rscript()
    if not exe:
        return None, "Rscript not found (set SCPCIP_RSCRIPT or use the scPCIP container image)"
    code = ("for (p in c(%s)) if (!requireNamespace(p, quietly = TRUE)) cat(p, '\\n')"
            % ", ".join(f"'{p}'" for p in R_PACKAGES))
    try:
        res = subprocess.run([exe, "--vanilla", "-e", code], capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, f"Rscript could not be run: {exc}"
    missing = res.stdout.split()
    if res.returncode != 0 or missing:
        return None, ("R package(s) not installed: " + ", ".join(missing)) if missing else res.stderr.strip()[-300:]
    return exe, ""


def organ_map(organ: str | None, species: str | None) -> str:
    if organ not in ("root", "leaf", "flower", "seed"):
        return "none"
    if organ == "leaf" and species and species.split()[0].capitalize() in GRASS_GENERA:
        return "leaf_grass"
    return organ


def _marker_genes(markers_path: Path) -> list[str]:
    import pandas as pd

    mk = pd.read_csv(markers_path, sep="\t", keep_default_na=False)
    genes: set[str] = set()
    for s in mk["query_member_genes"].astype(str):
        genes.update(g for g in s.split(",") if g)
    return sorted(genes)


def write_expression(adata: Any, cell_ids, genes: list[str], out_dir: Path) -> None:
    """log1p(counts per 10,000) of `genes` over the given cells (library size from all genes), genes x cells."""
    import scipy.io
    import scipy.sparse as sp

    names = np.asarray(adata.var_names.astype(str))
    pos = {g: i for i, g in enumerate(names)}
    keep = [g for g in genes if g in pos]
    idx = adata.obs_names.get_indexer(list(map(str, cell_ids)))
    if (idx < 0).any():
        raise ValueError("predictions contain cells that are not in the input")
    X = sp.csr_matrix(adata.X)[idx]
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    sub = X[:, [pos[g] for g in keep]].astype(np.float64)
    sub = sp.csr_matrix(sp.diags(1e4 / lib) @ sub)
    sub.data = np.log1p(sub.data)
    out_dir.mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(str(out_dir / "data.mtx"), sp.csc_matrix(sub.T))
    (out_dir / "genes.tsv").write_text("".join(f"{g}\n" for g in keep))
    (out_dir / "cells.tsv").write_text("".join(f"{c}\n" for c in map(str, cell_ids)))


def make(adata: Any, out_dir: str | Path, dataset_id: str, organ: str | None, species: str | None,
         exe: str | None = None) -> dict:
    """Draw the R figures for a finished predict run. `adata` = the gene-level query (raw counts). Returns
    {"status": "done", "dir": ..., "files": [...]} or {"status": "failed" | "skipped", "reason": ...}."""
    import pandas as pd

    out_dir = Path(out_dir)
    need = [out_dir / f"{dataset_id}.{k}.tsv" for k in ("predictions", "umap", "markers")]
    if not all(p.exists() for p in need):
        return {"status": "skipped", "reason": "run predict with the report first (needs " +
                ", ".join(p.name for p in need) + ")"}
    if exe is None:
        exe, why = check()
        if exe is None:
            return {"status": "skipped", "reason": why}
    fig_dir = out_dir / f"{dataset_id}.scPCIP_figures"
    cells = pd.read_csv(need[0], sep="\t", usecols=["cell_id"], dtype=str)["cell_id"]
    tmp = Path(tempfile.mkdtemp(prefix=".scpcip_r_", dir=out_dir))
    try:
        write_expression(adata, cells, _marker_genes(need[2]), tmp / "expr")
        for p in need:
            os.symlink(p.resolve(), tmp / p.name)
        omap = organ_map(organ, species)
        cmd = [exe, "--vanilla", str(script_path()), str(tmp), dataset_id, omap, str(fig_dir)]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
        except subprocess.TimeoutExpired:
            return {"status": "failed", "reason": f"R figures timed out after {TIMEOUT_S} s"}
        (out_dir / f"{dataset_id}.scPCIP_figures.log").write_text(res.stdout + res.stderr)
        if res.returncode != 0:
            tail = [ln for ln in res.stderr.strip().splitlines() if ln.strip()][-3:]
            return {"status": "failed", "reason": " | ".join(tail) or f"Rscript exit {res.returncode}",
                    "log": str(out_dir / f"{dataset_id}.scPCIP_figures.log")}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    files = sorted(str(p) for p in fig_dir.iterdir()) if fig_dir.exists() else []
    return {"status": "done", "dir": str(fig_dir), "organ_map": omap, "files": files}
