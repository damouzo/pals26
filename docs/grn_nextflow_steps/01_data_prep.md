# Step 1 — `CO_DATA_PREP` — Preprocessing & UMAP

**Process file:** `modules/celloracle_steps.nf`
**Script:** `bin/co_data_prep.py`

## What it does

Mirrors the first block of the CellOracle tutorial:

1. Downloads the scvelo pancreas dataset (`scv.datasets.pancreas()`)
2. Drops the pre-computed `highly_variable_genes` column and clears
   `uns / obsp / obsm`
3. Stores the raw counts in `adata.layers["counts"]` so CellOracle can
   re-use them later
4. Standard scanpy preprocessing:
   - `sc.pp.filter_genes(min_cells=1)`
   - `sc.pp.normalize_total` + `sc.pp.log1p`
   - `sc.pp.highly_variable_genes(subset=True)` to select HVGs
   - `sc.tl.pca(svd_solver='arpack')` → 50 PCs by default
   - `sc.pp.neighbors` and `sc.tl.leiden(resolution=0.3)`
   - `sc.tl.umap` for the 2-D embedding
5. Writes a Leiden-coloured UMAP to `qc_umap.png`
6. Pickles the AnnData to `adata_preprocessed.pkl` for the next
   process.

## Inputs

None — the dataset is fetched from scvelo's CDN on first use.

## Outputs

| File | Description |
|------|-------------|
| `adata_preprocessed.pkl` | Pickled AnnData with `counts` layer, HVGs, PCA, UMAP, Leiden |
| `qc_umap.png`            | Leiden-coloured UMAP (DPI 300) |
| `versions.yml`           | Tool versions for reproducibility |

## Knobs you may want to tune

| Flag                       | Default | Effect                       |
|----------------------------|---------|------------------------------|
| `--leiden-resolution`      | 0.3     | Cluster granularity          |
| `--filter-min-cells`       | 1       | Gene filter threshold        |

## Why pickle?

CellOracle's `Oracle.import_anndata_as_raw_count` needs every layer
and the full metadata. HDF5 is the alternative, but pickle round-trips
arbitrary Python objects (e.g. `uns` with custom dicts) without
serialisation boilerplate, which is ideal for a teaching pipeline.
