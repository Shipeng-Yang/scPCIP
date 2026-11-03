# Using scPCIP from an AI assistant

Two pieces, used together or separately:

| Piece | Path | Role |
|---|---|---|
| MCP server | `mcp/scpcip_mcp_server.py` | lets any MCP-capable assistant run scPCIP and read its results |
| Agent skill | `skills/scpcip/SKILL.md` | tells the assistant how to run scPCIP correctly and how to report the result |

## MCP server

The server runs the `scpcip` command line in your scPCIP environment, so it does not need scPCIP in its own
environment. Install the MCP SDK in a separate environment to keep the scPCIP environment unchanged:

```bash
python -m venv scpcip_mcp_venv
scpcip_mcp_venv/bin/pip install "mcp==2.3.0"
```

Register it with your assistant (stdio transport), for example in a JSON MCP configuration:

```json
{
  "mcpServers": {
    "scpcip": {
      "command": "/path/to/scpcip_mcp_venv/bin/python",
      "args": ["/path/to/scPCIP/integrations/mcp/scpcip_mcp_server.py"],
      "env": {"SCPCIP_PYTHON": "/path/to/scpcip_env/bin/python"}
    }
  }
}
```

Tools: `scpcip_supported`, `scpcip_inspect_input`, `scpcip_annotate`, `scpcip_markers`.

Check the installation end to end:

```bash
SCPCIP_PYTHON=/path/to/scpcip_env/bin/python \
  scpcip_mcp_venv/bin/python integrations/mcp/smoke_test.py data.h5ad "Zea mays" root ./scpcip_smoke
```

## Skill

Copy `skills/scpcip/` into the skills directory of your assistant. The skill works with the MCP server or, without
it, with the `scpcip` command line.
