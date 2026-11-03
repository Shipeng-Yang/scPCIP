"""Protocol test of the MCP server, independent of any assistant: every tool, normal and error inputs.

Usage: SCPCIP_PYTHON=... python protocol_test.py <gene-level root .h5ad> <species> <work dir>
Exits non-zero on the first failed check.
"""
import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

INP, SPECIES, WORK = sys.argv[1], sys.argv[2], Path(sys.argv[3])
RESULTS = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


async def call(s: ClientSession, tool: str, args: dict) -> dict:
    res = await s.call_tool(tool, args)
    return json.loads(res.content[0].text)


async def main() -> None:
    params = StdioServerParameters(command=sys.executable, args=[str(Path(__file__).with_name("scpcip_mcp_server.py"))],
                                   env=dict(os.environ, SCPCIP_OUT=str(WORK / "default_out")))
    async with stdio_client(params) as (r, w), ClientSession(r, w) as s:
        await s.initialize()
        names = sorted(t.name for t in (await s.list_tools()).tools)
        check("lists four tools", names == ["scpcip_annotate", "scpcip_inspect_input", "scpcip_markers",
                                            "scpcip_supported"], str(names))
        sup = await call(s, "scpcip_supported", {})
        check("supported organs", sup["organs"] == ["root", "leaf", "flower", "seed"])
        check("stem flagged unsupported", "stem" in sup["unsupported_organs"])
        check("species list non-empty", SPECIES in sup["species_with_shipped_gene_map"])
        ins = await call(s, "scpcip_inspect_input", {"input_path": INP})
        check("inspect reports shape", ins.get("n_cells", 0) > 0 and ins.get("n_genes", 0) > 0, str(ins.get("n_cells")))
        bad = await call(s, "scpcip_inspect_input", {"input_path": str(WORK / "missing.h5ad")})
        check("inspect missing file -> error", "error" in bad)
        e1 = await call(s, "scpcip_annotate", {"input_path": str(WORK / "missing.h5ad"), "species": SPECIES})
        check("annotate missing file -> error", "error" in e1 and "not found" in e1["error"])
        e2 = await call(s, "scpcip_annotate", {"input_path": INP, "species": SPECIES, "organ": "stem"})
        check("annotate stem -> refused", "error" in e2 and "stem" in e2["error"])
        e3 = await call(s, "scpcip_annotate", {"input_path": INP, "species": "Oryza sativa" if SPECIES != "Oryza sativa"
                                               else "Zea mays", "organ": "root", "out_dir": str(WORK / "wrong")})
        check("annotate wrong species -> failed run reported", "error" in e3 and "orthogroup" in e3.get("log_tail", ""))
        ok = await call(s, "scpcip_annotate", {"input_path": INP, "species": SPECIES, "organ": "root",
                                               "out_dir": str(WORK / "ok")})
        part = (ok.get("parts") or [{}])[0]
        check("annotate returns composition", bool(part.get("composition")), f"{part.get('n_cells')} cells")
        check("annotate writes HTML report", bool(part.get("html_report")) and Path(part["html_report"]).exists())
        mk = await call(s, "scpcip_markers", {"markers_table": part["markers_table"], "top": 2})
        check("markers per type", len(mk.get("markers", {})) >= 2 and all(len(v) <= 2 for v in mk["markers"].values()))
        one = await call(s, "scpcip_markers", {"markers_table": part["markers_table"], "cell_type": "endodermis"})
        check("markers single type filter", list(one.get("markers", {})) == ["endodermis"])
    failed = [n for n, ok, _ in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


asyncio.run(main())
