"""Locations of the shipped configs, dictionaries and the gene -> orthogroup map.

Configs and dictionaries are package data. The gene -> orthogroup map (about 50 MB) is downloaded once from the
Hugging Face Hub into a local cache; set SCPCIP_FEATURE_MAP to use a local copy instead.
"""

from __future__ import annotations

import os
from importlib import resources as _ir
from pathlib import Path

HF_REPO_ID = "shipengyang/scPCIP"
HF_REPO_TYPE = "model"
FEATURE_MAP_FILE = "gene_to_feature.tsv"
# Hub commit of the gene -> orthogroup map this release was validated with. A release that changes the map moves
# this pin, so an installed version always reads the same map.
FEATURE_MAP_REVISION = "7eaec1bb7081b0562f315e9c8389496807a5f137"
ORGANS = ("root", "leaf", "flower", "seed")


def data_dir() -> Path:
    return Path(str(_ir.files("pcip") / "data"))


def config_path(name: str) -> Path:
    return data_dir() / "configs" / name


def resolve(path: str | Path, base: Path) -> Path:
    """A path written in a config: absolute, relative to the config file, or relative to the package data."""
    p = Path(path)
    if p.is_absolute():
        return p
    for root in (base, data_dir()):
        if (root / p).exists():
            return root / p
    return base / p


def cache_dir() -> Path:
    root = os.environ.get("SCPCIP_CACHE") or os.path.join(os.environ.get("XDG_CACHE_HOME", "~/.cache"), "scpcip")
    return Path(root).expanduser()


def feature_map_path() -> Path:
    """Local path of the gene -> orthogroup map, downloading it on first use."""
    env = os.environ.get("SCPCIP_FEATURE_MAP")
    if env:
        return Path(env).expanduser()
    local = cache_dir() / FEATURE_MAP_REVISION[:12] / FEATURE_MAP_FILE
    if local.exists():
        return local
    try:
        from huggingface_hub import hf_hub_download
    except ImportError as exc:  # pragma: no cover - dependency is declared
        raise SystemExit("scpcip: huggingface_hub is needed to download the gene -> orthogroup map "
                         "(pip install huggingface_hub), or set SCPCIP_FEATURE_MAP to a local copy") from exc
    print(f"scpcip: downloading the gene -> orthogroup map from huggingface.co/{HF_REPO_ID} (once) ...")
    path = hf_hub_download(repo_id=HF_REPO_ID, repo_type=HF_REPO_TYPE, filename=FEATURE_MAP_FILE,
                           revision=FEATURE_MAP_REVISION, local_dir=str(local.parent))
    return Path(path)


def dictionary_names() -> dict[str, str]:
    """organ -> file name of the dictionary its shipped config uses."""
    import yaml

    names = {}
    for organ in ORGANS:
        with config_path(f"predict_{organ}.yaml").open(encoding="utf-8") as fh:
            names[organ] = Path(yaml.safe_load(fh)["dictionary"]).stem
    return names


def version_text(version: str) -> str:
    dicts = ", ".join(dictionary_names())
    return (f"scpcip {version}\n  dictionaries: {dicts}\n"
            f"  gene -> orthogroup map: huggingface.co/{HF_REPO_ID} @ {FEATURE_MAP_REVISION[:12]}")
