#!/usr/bin/env python
"""
co_perturbation.py — Step 4 of the pals26 pipeline.

Takes the trained Oracle + Links objects, perturbs a transcription
factor (default: knock-out of Mafb), and propagates the perturbation
through the inferred GRN. Saves the perturbed Oracle for the final
gradient / development module step.
"""
from __future__ import annotations

import argparse
import pickle
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path
import celloracle as co


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-adata",      required=True)
    p.add_argument("--input-oracle",     required=True)
    p.add_argument("--input-links",      required=True)
    p.add_argument("--cluster-col",      default="leiden")
    p.add_argument("--perturb-tf",       default="Mafb")
    p.add_argument("--perturb-value",    type=float, default=0.0)
    p.add_argument("--n-propagation",    type=int,   default=3)
    p.add_argument("--n-neighbors",      type=int,   default=200)
    p.add_argument("--sigma-corr",       type=float, default=0.05)
    p.add_argument("--p-mass-n-grid",    type=int,   default=40)
    p.add_argument("--p-mass-smooth",    type=float, default=0.8)
    p.add_argument("--p-mass-min-mass",  type=float, default=6.2)
    p.add_argument("--grn-alpha",        type=float, default=10.0)
    p.add_argument("--output-oracle",    required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    with open(args.input_oracle, "rb") as fh:
        oracle: co.Oracle = pickle.load(fh)
    with open(args.input_links, "rb") as fh:
        links = pickle.load(fh)

    # Cluster-specific TF dictionary (only TFs with edges in the cluster)
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)

    # Fit a regression model per cluster for in-silico simulation
    # alpha is the ridge regularisation strength (CellOracle uses the
    # same `alpha` symbol as the GRN fit; default 10 follows the tutorial).
    oracle.fit_GRN_for_simulation(
        alpha=args.grn_alpha,
        use_cluster_specific_TFdict=True,
    )

    # Simulate the actual perturbation
    oracle.simulate_shift(
        perturb_condition={args.perturb_tf: args.perturb_value},
        n_propagation=args.n_propagation,
    )

    # Estimate transition probabilities & embedding shift
    oracle.estimate_transition_prob(
        n_neighbors=args.n_neighbors,
        knn_random=True,
        sampled_fraction=1,
    )
    oracle.calculate_embedding_shift(sigma_corr=args.sigma_corr)

    # Build the simulation grid + mass filter
    oracle.calculate_p_mass(
        smooth=args.p_mass_smooth,
        n_grid=args.p_mass_n_grid,
        n_neighbors=args.n_neighbors,
    )
    oracle.calculate_mass_filter(
        min_mass=args.p_mass_min_mass,
        plot=False,
    )

    #Plots
    # 1. UMAP plot of imputed Mafb expression
    sc.pl.umap(oracle.adata, color=[args.perturb_tf], layer="imputed_count", show=False)
    plt.savefig(Path(args.output_oracle).parent / f"umap_imputed_{args.perturb_tf}.png", dpi=300)
    plt.close()

    # 2. Real vs randomized simulation arrows (Quiver Plots)
    fig, ax = plt.subplots(1, 2, figsize=[13, 6])
    scale = 30 
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {args.perturb_tf} KO")
    
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title("Randomized simulation vector")
    
    plt.tight_layout()
    plt.savefig(Path(args.output_oracle).parent / f"perturbation_quiver_{args.perturb_tf}.png", dpi=300)
    plt.close()

    # Re-serialise the perturbed Oracle for the gradient step
    with open(args.output_oracle, "wb") as fh:
        pickle.dump(oracle, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[co_perturbation] simulated {args.perturb_tf} = {args.perturb_value} "
          f"({args.n_propagation} propagation steps)")


if __name__ == "__main__":
    main()
