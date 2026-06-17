#!/usr/bin/env python
"""
co_pseudotime.py — Step 2 of the pals26 pipeline.

Loads the pre-processed AnnData, infers a diffusion-pseudotime
trajectory using CellOracle's Pseudotime_calculator, and writes
back the AnnData with the `Pseudotime` obs column attached.

Cluster column name, lineage name, and root-cell selection are all
parameterised so the script can be reused on other datasets.
"""
from __future__ import annotations

import argparse
import pickle

import numpy as np

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from celloracle.applications import Pseudotime_calculator


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-adata",   required=True)
    p.add_argument("--cluster-col",   default="leiden",
                   help="obs column carrying the cluster labels")
    p.add_argument("--lineage-name",  default="the_one",
                   help="name of the lineage that groups all clusters")
    p.add_argument("--root-strategy", default="max-umap-sum",
                   choices=["max-umap-sum", "first-cluster"],
                   help="How to pick the root cell automatically")
    p.add_argument("--root-cluster",  default=None,
                   help="When root-strategy=first-cluster, name of the cluster to use")
    p.add_argument("--output-adata",  required=True)
    p.add_argument("--output-plot",   required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    with open(args.input_adata, "rb") as fh:
        adata: ad.AnnData = pickle.load(fh)

    # CellOracle's Pseudotime_calculator operates on a copy of adata
    pt = Pseudotime_calculator(
        adata=adata,
        obsm_key="X_umap",
        cluster_column_name=args.cluster_col,
    )

    # One lineage containing every cluster, mirroring the tutorial
    lineage_dict = {args.lineage_name: adata.obs[args.cluster_col].unique().tolist()}
    pt.set_lineage(lineage_dictionary=lineage_dict)

    # Root-cell selection
    if args.root_strategy == "max-umap-sum":
        coord_sum = np.asarray(adata.obsm["X_umap"]).sum(axis=1)
        max_idx = int(np.argmax(coord_sum))
        root_cells = {args.lineage_name: adata.obs_names[max_idx]}
    else:
        if not args.root_cluster:
            raise SystemExit("ERROR: --root-cluster is required when --root-strategy=first-cluster")
        mask = adata.obs[args.cluster_col] == args.root_cluster
        if not mask.any():
            raise SystemExit(
                f"ERROR: cluster '{args.root_cluster}' not present in obs['{args.cluster_col}']"
            )
        root_cells = {args.lineage_name: adata.obs_names[mask.argmax()]}
    pt.set_root_cells(root_cells=root_cells)

    # Diffusion pseudotime — needs a graph on the UMAP / PCA space
    sc.tl.diffmap(pt.adata)
    pt.get_pseudotime_per_each_lineage()

    # Plot & persist
    pt.plot_pseudotime(cmap="rainbow")
    ax = plt.gca()
    if ax.collections:
        scatter = ax.collections[0]
        plt.colorbar(scatter, ax=ax, label="Pseudotime")
    plt.tight_layout()
    plt.savefig(args.output_plot, dpi=300)
    plt.close()

    # Carry the pseudotime back to the original adata and re-serialise
    adata.obs["Pseudotime"] = pt.adata.obs["Pseudotime"]
    with open(args.output_adata, "wb") as fh:
        pickle.dump(adata, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[co_pseudotime] pseudotime range: "
          f"{adata.obs['Pseudotime'].min():.3f} -> {adata.obs['Pseudotime'].max():.3f}")


if __name__ == "__main__":
    main()
