"""End-to-end check: start the server over stdio, list its tools and annotate one dataset.

Usage: python smoke_test.py <input.h5ad> <species> <organ> <out_dir>
"""
import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client


async def main(inp: str, species: str, organ: str, out: str) -> None:
    params = StdioServerParameters(command=sys.executable, args=[str(Path(__file__).with_name("scpcip_mcp_server.py"))],
                                   env=dict(os.environ))
    async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
        await session.initialize()
        tools = await session.list_tools()
        print("tools:", [t.name for t in tools.tools])
        for name, args in [("scpcip_supported", {}), ("scpcip_inspect_input", {"input_path": inp}),
                           ("scpcip_annotate", {"input_path": inp, "species": species, "organ": organ,
                                                "out_dir": out})]:
            res = await session.call_tool(name, args)
            text = res.content[0].text
            print(f"--- {name} ---\n{text[:3000]}")
            if name == "scpcip_annotate":
                summary = json.loads(text)
                mk = summary["parts"][0]["markers_table"]
                res = await session.call_tool("scpcip_markers", {"markers_table": mk, "top": 3})
                print(f"--- scpcip_markers ---\n{res.content[0].text[:2000]}")


asyncio.run(main(*sys.argv[1:5]))
