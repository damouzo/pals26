# GRN with your own scMultiome

A short guide to inferring a Gene Regulatory Network with **CellOracle** when you have paired **scRNA-seq + scATAC-seq** from the same cells. You start from a Seurat object that already has both assays and a `CellType` annotation; three scripts take it from there to a per-cluster GRN.

## What you need

- A Seurat object with `RNA` and `ATAC` assays built from the same cells
- A `CellType` column in `meta.data`.
- A `parameters_celloracle.yaml` with these keys:
  - `Path.Home` — project root.
  - `Analysis_Name` — subfolder name; raw data and outputs will live under it.
  - `Path.SeuObj` — path to the input Seurat `.rds`.
  - `Info4GRN.clustering` — column with cell-type labels.
  - `Info4GRN.umap` — UMAP embedding key (e.g. `umap_multiome`).
  - `Info4GRN.clust2study` — single cluster to export as its own CSV.
  - `Path.genome_ref` — CellOracle genome reference directory.

## The flow

1. **Integrate** — `Multiome_Integration.R` joins RNA and ATAC into one weighted graph and writes a new Seurat object.
2. **Prepare inputs** — `GRN_dataProcess.R` exports RNA and ATAC as `.h5ad`, runs **Cicero** for peak coaccessibility.
3. **Build the GRN** — `GRN_analysis.py` runs CellOracle: motif scan on coaccessible peaks, KNN imputation, per-cluster ridge regression, network scores.

## Step 1 — Integrate

Run `Multiome_Integration.R`. It expects a `Parameters` file (key:value lines) and a Seurat object already annotated with `CellType` and Harmony embeddings for both modalities.

It performs `FindMultiModalNeighbors` over `harmony_rna` + `harmony_atac`, runs UMAP on the weighted graph (`umap_multiome`), and re-clusters at resolution 0.5.

Output: `SeuObj_MultiomeIntegr.rds`.

## Step 2 — Prepare inputs

Run `GRN_dataProcess.R` with the same Seurat object. It does three things:

- Exports the RNA assay to `HSC_5BM_RNA.h5ad` and the ATAC assay to `HSC_5BM_ATAC.h5ad` (see note in *Notes* about renaming these).
- Converts the ATAC assay to a `cell_data_set`, preprocesses with LSI, and runs **Cicero** to get peak–peak coaccessibility.
- Writes `cicero_connections.csv`, `all_peaks.csv`, and two PDF histograms (`Histogram_Cicero.pdf`, `Histogram_Coaccess.pdf`) for sanity checks.

Look at the two histograms before continuing. A bimodal or skewed coaccess distribution is normal; flat or empty means Cicero did not converge.

## Step 3 — Build the GRN

Run `GRN_analysis.py` with `parameters_celloracle.yaml` in the working directory. The script:

- Loads the ATAC `.h5ad` and Cicero connections; annotates TSSs and keeps peaks with `coaccess >= 1`.
- Scans TF binding motifs with the `gimme.vertebrate.v5.0` database to build the base GRN.
- Loads the RNA `.h5ad`, builds a CellOracle `Oracle`, runs PCA and KNN imputation (`k ≈ 0.025 × n_cells`).
- Fits a ridge regression per cluster (`alpha = 10`) and exports one CSV per cluster plus an Excel with all clusters.
- Filters edges (`p = 0.001`, top 2000 by `|coef|`) and computes network scores (degree centrality, in/out/all).

Outputs:
- `base_GRN_dataframe_gimme.vertebrate.v5.0.parquet`
- `raw_gimme.vertebrate.v5.0_GRN_for_<clust2study>_.csv`
- `raw_gimme.vertebrate.v5.0_GRN_clusters.xlsx`
- `links_merged_score.xlsx`
- Plots under `output/plots/`.

The motif scan and the per-cluster ridge fit are the slowest steps. Expect tens of minutes on a laptop for a few thousand cells.

