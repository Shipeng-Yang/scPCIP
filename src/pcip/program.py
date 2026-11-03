"""Program dictionary loading.

A program dictionary is a YAML description of weighted gene programs over orthogroups (one or more programs per
cell type, plus optional fine-type and developmental-state programs).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_dictionary(path: str | Path) -> dict[str, Any]:
    """Load a program dictionary YAML file into a plain dict."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"dictionary not found: {path}")
    with path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"dictionary root must be a mapping, got {type(data).__name__}")
    return data
