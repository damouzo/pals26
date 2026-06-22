#!/usr/bin/env python
"""
co_gradient_dev.py — Step 5 of the pals26 pipeline.

Loads the perturbed Oracle, attaches the pseudotime computed in step 2,
builds a Gradient_calculator on the same grid as the perturbation
simulation, runs Oracle_development_module, and saves the final
inner-product-on-grid plot + a per-cell dev-score table.
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd

import celloracle as co
from celloracle.applications import Gradient_calculator, Oracle_development_module


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-adata",    required=True)
    p.add_argument("--input-oracle",   required=True)
    p.add_argument("--pseudotime-key", default="Pseudotime")
    p.add_argument("--n-grid",         type=int,   default=40)
    p.add_argument("--n-neighbors",    type=int,   default=200)
    p.add_argument("--smooth",         type=float, default=0.8)
    p.add_argument("--min-mass",       type=float, default=6.2)
    p.add_argument("--n-poly",         type=int,   default=3)
    p.add_argument("--n-bins",         type=int,   default=10)
    p.add_argument("--output-plot",    required=True)
    p.add_argument("--output-scores",  required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Load perturbed Oracle and inject the pseudotime we computed earlier
    with open(args.input_oracle, "rb") as fh:
        oracle: co.Oracle = pickle.load(fh)
    with open(args.input_adata, "rb") as fh:
        adata_pt: ad.AnnData = pickle.load(fh)
    oracle.adata.obs[args.pseudotime_key] = adata_pt.obs[args.pseudotime_key]

    # --- Gradient calculator (identical grid to the perturbation) ---
    gradient = Gradient_calculator(
        oracle_object=oracle,
        pseudotime_key=args.pseudotime_key,
    )
    gradient.calculate_p_mass(
        smooth=args.smooth,
        n_grid=args.n_grid,
        n_neighbors=args.n_neighbors,
    )
    gradient.calculate_mass_filter(min_mass=args.min_mass, plot=False)
    gradient.transfer_data_into_grid(
        args={"method": "polynomial", "n_poly": args.n_poly},
        plot=False,
    )
    gradient.calculate_gradient()

    # --- Development module: perturbation × reference flow ---
    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=args.n_bins)

    # --- Final perturbation score plot ---
    fig, ax = plt.subplots(figsize=[6, 6])
    dev.plot_inner_product_on_grid(vm=1, s=50, ax=ax)
    dev.plot_simulation_flow_on_grid(scale=30, show_background=False, ax=ax)
    norm = colors.Normalize(vmin=-1, vmax=1)
    sm = cm.ScalarMappable(cmap="PiYG_r", norm=norm) # PiYG_r: green (+), pink (-)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Inner Product Score (Green: Pro-development | Pink: Anti-development )")
    plt.tight_layout()
    plt.savefig(args.output_plot, dpi=300)
    plt.close()

    # --- Per-cell digitised inner-product scores (TSV) ---
    scores = dev.digitized_ip_df if hasattr(dev, "digitized_ip_df") else pd.DataFrame()
    scores.to_csv(args.output_scores, sep="\t")

    print(f"[co_gradient_dev] wrote {args.output_plot} and {args.output_scores}")


if __name__ == "__main__":
    main()
