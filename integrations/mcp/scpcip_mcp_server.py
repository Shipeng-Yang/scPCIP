"""MCP server exposing scPCIP cell-type annotation to MCP-capable assistants.

The server does not import scPCIP. It runs the `scpcip` command line in the environment where scPCIP is installed
and returns compact summaries, so the assistant never has to read large tables.

Environment variables
  SCPCIP_PYTHON   Python of the environment where scPCIP is installed (default: python3 on PATH)
  SCPCIP_OUT      default parent directory for results (default: ./scpcip_runs)

Run: python scpcip_mcp_server.py            (stdio transport)
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import time
from pathlib import Path

from mcp.server import MCPServer

PY = os.environ.get("SCPCIP_PYTHON", "python3")
OUT_ROOT = Path(os.environ.get("SCPCIP_OUT", "scpcip_runs")).resolve()
ORGANS = ("root", "leaf", "flower", "seed")

server = MCPServer(
    name="scpcip",
    instructions=(
        "scPCIP annotates plant single-cell data (root, leaf, flower, seed) with cross-species gene programs. "
        "Ask the user for species and organ if they are not stated; never guess the organ from expression. "
        "Every cell receives a type (the tool never abstains); report the confidence quartiles with the result."
    ),
)


def _run(args: list[str], timeout: int | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([PY, *args], capture_output=True, text=True, timeout=timeout)


def _read_tsv(path: Path, limit: int | None = None) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[:limit] if limit else rows


def _species_in_map() -> list[str]:
    res = _run(["-m", "pcip", "species"], timeout=1800)
    if res.returncode != 0:
        return []
    return [ln.split("\t", 1)[0] for ln in res.stdout.splitlines() if "\t" in ln]


@server.tool()
def scpcip_supported() -> dict:
    """List the organs and species scPCIP can annotate out of the box."""
    return {
        "organs": list(ORGANS),
        "unsupported_organs": {"stem": "not shipped (below the single-marker baseline on held-out species)"},
        "species_with_shipped_gene_map": _species_in_map(),
        "other_species": "supply your own gene-to-orthogroup table (eggNOG-mapper output) as feature_map",
        "input": "AnnData .h5ad with raw counts, genes as var_names (gene IDs of the species' reference annotation)",
    }


@server.tool()
def scpcip_inspect_input(input_path: str) -> dict:
    """Describe an .h5ad file: cell and gene counts, example gene IDs and the cell metadata columns."""
    code = (
        "import anndata as ad, json, sys\n"
        "a = ad.read_h5ad(sys.argv[1], backed='r')\n"
        "obs = {c: [str(v) for v in a.obs[c].astype(str).unique()[:8]] for c in a.obs.columns[:40]}\n"
        "print(json.dumps({'n_cells': a.n_obs, 'n_genes': a.n_vars, 'example_genes': list(map(str, a.var_names[:8])),"
        " 'obs_columns': obs}))\n"
    )
    res = _run(["-c", code, str(Path(input_path).expanduser().resolve())], timeout=600)
    if res.returncode != 0:
        return {"error": res.stderr.strip()[-2000:]}
    return json.loads(res.stdout.strip().splitlines()[-1])


@server.tool()
def scpcip_annotate(input_path: str, species: str, organ: str = "root", organ_key: str | None = None,
                    organ_map: str | None = None, feature_map: str | None = None, og_level: bool = False,
                    out_dir: str | None = None, report: bool = True) -> dict:
    """Annotate every cell of an .h5ad file with one cell type and its tissue system.

    organ: root, leaf, flower or seed for a single-organ file. For a file mixing organs, give organ_key (the cell
    metadata column holding the organ, or the sample when organ_map maps samples to organs) instead.
    feature_map: gene-to-orthogroup TSV for a species without a shipped map.
    report: also write figures, marker tables and a self-contained HTML report.
    """
    src = Path(input_path).expanduser().resolve()
    if not src.exists():
        return {"error": f"input not found: {src}"}
    if organ_key is None and organ not in ORGANS:
        return {"error": f"organ must be one of {', '.join(ORGANS)} (stem is not supported)"}
    out = Path(out_dir).expanduser().resolve() if out_dir else OUT_ROOT / f"{src.stem}_{time.strftime('%Y%m%d_%H%M%S')}"
    args = ["-m", "pcip", "predict", "--input", str(src), "--species", species, "--out", str(out)]
    if organ_key:
        args += ["--organ-key", organ_key] + (["--organ-map", organ_map] if organ_map else [])
    else:
        args += ["--organ", organ]
    if feature_map:
        args += ["--feature-map", str(Path(feature_map).expanduser().resolve())]
    if og_level:
        args.append("--og-level")
    args.append("--report" if report else "--no-report")
    t0 = time.time()
    res = _run(args)
    log = (res.stdout + res.stderr).strip()
    if res.returncode != 0:
        return {"error": "scpcip predict failed", "exit_code": res.returncode, "log_tail": log[-3000:]}
    return {"out_dir": str(out), "seconds": round(time.time() - t0, 1),
            "warnings": [ln for ln in log.splitlines() if "WARNING" in ln], **_summarise(out)}


def _summarise(out: Path) -> dict:
    reports = sorted(out.rglob("*.predict_report.json"))
    if not reports:
        return {"note": "no predict_report.json found", "files": sorted(p.name for p in out.iterdir())}
    parts = []
    for rp in reports:
        meta = json.loads(rp.read_text())
        stem = rp.name[: -len(".predict_report.json")]
        comp = rp.with_name(f"{stem}.reliability.tsv")
        comp = comp if comp.exists() else rp.with_name(f"{stem}.composition.tsv")
        parts.append({
            "dataset_id": meta.get("dataset_id"), "n_cells": meta.get("n_cells"),
            "dictionary": Path(str(meta.get("dictionary", ""))).name,
            "composition": [{k: _num(v) for k, v in r.items()} for r in _read_tsv(comp)] if comp.exists() else None,
            "confidence_note": meta.get("confidence"),
            "html_report": str(rp.with_name(f"{stem}.report.html")) if rp.with_name(f"{stem}.report.html").exists() else None,
            "predictions_table": str(rp.with_name(f"{stem}.predictions.tsv")),
            "markers_table": str(rp.with_name(f"{stem}.markers.tsv")) if rp.with_name(f"{stem}.markers.tsv").exists() else None,
        })
    return {"parts": parts}


def _num(v: str):
    try:
        f = float(v)
        return int(f) if f.is_integer() and re.fullmatch(r"-?\d+", v) else round(f, 4)
    except (TypeError, ValueError):
        return v


@server.tool()
def scpcip_markers(markers_table: str, cell_type: str | None = None, top: int = 10) -> dict:
    """Top-ranked marker genes per predicted cell type from a markers table written by scpcip_annotate."""
    rows = _read_tsv(Path(markers_table).expanduser())
    if not rows:
        return {"error": "empty markers table"}
    key = next((k for k in ("cell_type", "predicted", "group") if k in rows[0]), list(rows[0])[0])
    out: dict[str, list] = {}
    for r in rows:
        ct = r[key]
        if cell_type and ct != cell_type:
            continue
        if len(out.setdefault(ct, [])) < top:
            out[ct].append({k: _num(v) for k, v in r.items() if k != key})
    return {"cell_type_column": key, "markers": out}


if __name__ == "__main__":
    server.run("stdio")
