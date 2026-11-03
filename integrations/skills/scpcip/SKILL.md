---
name: scpcip
description: Annotate plant single-cell RNA-seq data (root, leaf, flower or seed; any species with a gene-to-orthogroup map) with scPCIP cross-species gene programs, and explain the result. Use when a user asks to annotate, label or identify cell types or tissues in a plant .h5ad file, or asks for marker genes of plant cell types from such data.
---

# scPCIP annotation

scPCIP assigns every cell one cell type and its tissue system by scoring cross-species gene programs (weighted
eggNOG orthogroups learned from many labelled plant atlases). It needs no reference atlas from the user's species.

## Before running

1. **Species and organ.** Both are required. If the user has not stated them, ask. Never infer the organ from
   expression. Supported organs: root, leaf, flower, seed. Stem is not supported; say so and stop.
2. **Mixed organs.** If one file holds several organs, ask which cell metadata column records the organ (or the
   sample, plus a sample-to-organ table) and use `organ_key` / `organ_map` instead of `organ`.
3. **Input.** An AnnData `.h5ad` with raw counts and gene IDs of the species' reference annotation as `var_names`.
   Run `scpcip_inspect_input` first: check the gene IDs look like the species' IDs and note the metadata columns.
4. **Species without a shipped map.** `scpcip_supported` lists species with a shipped gene-to-orthogroup map. For
   another species the user must supply one (eggNOG-mapper output converted to the scPCIP table) as `feature_map`.

## Running

Call `scpcip_annotate` with the input path, species and organ. Keep `report` on: it writes figures, tables and a
self-contained HTML report. Without the MCP server, run the same thing on the command line with absolute paths:

```bash
scpcip predict --input /abs/data.h5ad --species "Oryza sativa" --organ root --out /abs/my_run
# organs: root, leaf, flower, seed
# mixed organs: replace --organ with --organ-key <column>  (or --organ-key sample --organ-map sample_to_organ.tsv)
# species with a shipped gene map: scpcip species
```

If the run stops because no gene maps to an orthogroup, the gene IDs do not belong to the stated species or its
annotation version. Tell the user what the IDs look like; do not report any annotation from that run.

## Reporting the result

- Give the composition table (cell type, cells, fraction) and the path to the HTML report.
- **Every cell gets a type; the tool never abstains.** Always state the confidence: `confidence_quartile` splits cells
  by margin within the query, and `expected_accuracy` is the measured accuracy of that quartile on held-out species.
  Treat quartile-1 cells as tentative.
- `cross_species_f1` in the composition is each type's measured accuracy on unseen species; flag types with low
  values rather than presenting all calls as equally reliable.
- If the log warns of a low gene-to-orthogroup mapping rate, say so: the gene IDs may not match the map.
- A cell type absent from the sample can still be assigned (mostly in quartile 1); mention this when a type is
  biologically unexpected.
- Only coarse types are reported. Do not present them as fine subtypes.
- Use `scpcip_markers` for marker genes: each row is a program orthogroup with its member genes in the user's
  species and their enrichment in the predicted type. Weight-ranked orthogroups are program evidence, not proof of
  expression in that dataset; prefer rows with high `detect_frac_in_type` and `log2fc_in_vs_rest`.

## Do not

- Do not rerun with a different organ to get a "better" answer.
- Do not describe cell-cycle, ribosomal or stress signals as cell identity.
- Do not claim accuracy beyond the numbers in the output files.
