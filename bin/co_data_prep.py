#!/usr/bin/env python
"""
co_data_prep.py — Step 1 of the pals_celloracle pipeline.

Loads the scvelo pancreas dataset, normalises, log-transforms, selects
highly-variable genes, computes PCA / neighbours / Leiden / UMAP, and
saves the resulting AnnData as a pickle for downstream processes.

"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use("Agg")  # headless rendering inside the container
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset",            default="pancreas",
                   help="scv.datasets.<name>() entry point (default: pancreas)")
    p.add_argument("--leiden-resolution",  type=float, default=0.3)
    p.add_argument("--filter-min-cells",   type=int,   default=1)
    p.add_argument("--output-adata",       required=True)
    p.add_argument("--output-umap",        required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # 1. Load the dataset dynamically from scvelo
    loader = getattr(scv.datasets, args.dataset)
    adata: ad.AnnData = loader()

    # 2. Strip pre-computed artefacts to keep the object minimal
    if "highly_variable_genes" in adata.var.columns:
        adata.var = adata.var.drop(columns=["highly_variable_genes"])
    adata.uns, adata.obsp, adata.obsm = {}, {}, {}

    # 3. Preserve raw counts so CellOracle can re-derive them later
    adata.layers["counts"] = adata.X.copy()

    # 4. Standard scanpy preprocessing
    sc.pp.filter_genes(adata, min_cells=args.filter_min_cells)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, subset=True, inplace=True)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=args.leiden_resolution)
    sc.tl.umap(adata)

    # 5. QC plot (clustered UMAP)
    sc.pl.umap(adata, color="leiden", legend_loc="on data",
               show=False, frameon=False)
    plt.tight_layout()
    plt.savefig(args.output_umap, dpi=300)
    plt.close()

        #celltype humap
    if "clusters_coarse" in adata.obs.columns:
        sc.pl.umap(adata, color="clusters_coarse", legend_loc="on data",
                   show=False, frameon=False)
        plt.tight_layout()
        plt.savefig(args.output_umap.replace("qc_umap.png", "celltype_umap.png"), dpi=300)
        plt.close()

    # 6. Serialise — pickle preserves every layer / obsm / uns the
    #    later processes need (Oracle.adata, Links, etc.)
    with open(args.output_adata, "wb") as fh:
        pickle.dump(adata, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[co_data_prep] saved {adata.n_vars} genes x {adata.n_obs} cells "
          f"({adata.obs['leiden'].nunique()} clusters) -> {args.output_adata}")


if __name__ == "__main__":
    main()
