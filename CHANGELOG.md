# Changelog

Version numbers follow `MAJOR.MINOR.PATCH`; see "Versions and updates" in the README for what each kind of release
may change.

## [1.0.0] - 2026-11-03

First release.

- Cell-type and tissue-system annotation of root, leaf, flower and seed single-cell RNA-seq data from any plant
  species with a gene -> orthogroup map; one program dictionary per organ.
- Per-cell program scores, margin-based confidence quartile and expected accuracy measured on held-out species.
- Figures, marker tables and a self-contained HTML report; optional R figures (organ map, UMAPs, marker heatmap and
  violins) and `scpcip plot`.
- Reproducibility check for multi-sample data with integrated clusters.
- `scpcip build-map` for species not in the shipped map; Python interface `pcip.annotate`.
- MCP server and agent skill in `integrations/`; Apptainer / Singularity image definition in `container/`.
