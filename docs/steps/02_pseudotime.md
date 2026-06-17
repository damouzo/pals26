# Step 2 — `CO_PSEUDOTIME` — Diffusion pseudotime

**Process file:** `modules/celloracle_steps.nf`
**Script:** `bin/co_pseudotime.py`

## What it does

CellOracle's `Pseudotime_calculator` wraps a diffusion-map approach
to assign a pseudotime value to every cell. The pipeline:

1. Loads the pickled AnnData
2. Builds a `Pseudotime_calculator` on the UMAP embedding using the
   Leiden cluster labels
3. Defines a single lineage (`the_one`) that contains every Leiden
   cluster
4. Picks the cell with the largest sum of UMAP coordinates as the
   root (matches the upstream tutorial's heuristic)
5. Computes the diffusion map on the calculator's AnnData
   (`sc.tl.diffmap`)
6. Runs `get_pseudotime_per_each_lineage()`
7. Renders a rainbow-coloured pseudotime plot (`pt.png`)
8. Copies the resulting `Pseudotime` column back to the original
   AnnData and pickles it.

## Inputs

| Channel | File | Description |
|---------|------|-------------|
| `input-adata` | `adata_preprocessed.pkl` | Output of step 1 |

## Outputs

| File | Description |
|------|-------------|
| `adata_with_pseudotime.pkl` | AnnData with `obs['Pseudotime']` |
| `pt.png`                    | Pseudotime-coloured UMAP       |
| `versions.yml`              | Tool versions                  |

## Knobs

| Flag               | Default  | Effect                            |
|--------------------|----------|-----------------------------------|
| `--cluster-col`    | `leiden` | Cluster column to use for lineage |
| `--lineage-name`   | `the_one` | Name of the single lineage       |

## Why pick a root cell automatically?

The upstream tutorial chooses the cell with the maximum sum of UMAP
coordinates; this is a quick way to anchor the trajectory on a
plausible "progenitor" cluster for the pancreas dataset. For real
analyses you would pass `set_root_cells` with a biologically-motivated
list of cells.
