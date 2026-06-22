# Step 3 — `CO_GRN_INFERENCE` — GRN model

**Process file:** `modules/celloracle_steps.nf`
**Script:** `bin/co_grn_inference.py`

## What it does

This is the heart of the pipeline. It builds a `co.Oracle` object,
imports the mouse scATAC-atlas base GRN, performs KNN-imputation of
dropouts, and fits a per-cluster ridge regression against every
transcription factor in the base GRN.

1. Loads the AnnData and restores the **raw count** layer
   (`adata.X = adata.layers['counts']`) — CellOracle needs un-logged
   counts.
2. Downloads the mouse scATAC base GRN via
   `co.data.load_mouse_scATAC_atlas_base_GRN()` (cached in the Docker
   image at build time, so it does **not** require network access on
   Apocrita).
3. `Oracle.import_anndata_as_raw_count(cluster_column_name='leiden',
   embedding_name='X_umap')` — registers the AnnData with the Oracle.
4. `Oracle.import_TF_data(TF_info_matrix=base_GRN)` — attaches the
   candidate TFs and their candidate binding sites.
5. `Oracle.perform_PCA()` — recomputes PCA inside the Oracle.
6. `Oracle.knn_imputation(n_pca_dims=26, k=92, balanced=True,
   b_sight=k*8, b_maxl=k*4, n_jobs=4)` — fills in dropouts. The
   parameters (`k`, `b_sight`, `b_maxl`) come from the upstream
   tutorial, which targets ~0.025 × n_cells.
7. `Oracle.get_links(cluster_name_for_GRN_unit='leiden', alpha=10)`
   — fits a separate ridge model per cluster, yielding a `Links`
   object with one edge-list per cluster.
8. `Links.filter_links(p=0.001, weight='coef_abs',
   threshold_number=10000)` — keeps the 10 000 most-confident edges
   per cluster (cap chosen to match the tutorial).
9. `Links.get_network_score()` — computes the GRN-boosting score used
   by the perturbation module.
10. Pickles the Oracle and Links objects.

## Inputs

| Channel | File | Description |
|---------|------|-------------|
| `input-adata` | `adata_with_pseudotime.pkl` | Output of step 2 (pseudotime is not used here, but carrying it forward avoids re-loading) |

## Outputs

| File | Description |
|------|-------------|
| `oracle_object.pkl` | Trained Oracle (PCA, KNN-imputed expression) |
| `links_object.pkl`  | Per-cluster GRN edges and network scores  |
| `grn_score.png`     | Diagnostic plot of cluster-level scores  |
| `versions.yml`      | Tool versions                             |

## Knobs

| Flag                | Default | Effect                              |
|---------------------|---------|-------------------------------------|
| `--n-pca`           | 26      | PCA components used for imputation  |
| `--knn-k`           | 92      | k for KNN imputation                |
| `--knn-n-jobs`      | 4       | Threads for KNN search              |
| `--grn-alpha`       | 10      | Ridge α for GRN regression          |
| `--filter-p`        | 0.001   | p-value cutoff for `filter_links`   |
| `--filter-max`      | 10000   | Max edges kept per cluster          |

## Why is this the slowest step?

The per-cluster ridge regression scans every (TF × gene) pair in the
base GRN, which is ≈ a few hundred thousand pairs across the ~3 696
pancreas genes. On a workstation this takes 5-15 minutes depending on
the number of clusters.
