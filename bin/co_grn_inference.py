#!/usr/bin/env python
"""
co_grn_inference.py — Step 3 of the pals_celloracle pipeline.

Builds a CellOracle Oracle object from the pre-processed AnnData,
imports the mouse scATAC base GRN, performs KNN-imputation, fits
cluster-wise GRNs and writes the trained Oracle + Links objects to
disk for the perturbation step.

PCA dimensions and KNN k can be supplied explicitly *or* set to
``"auto"`` so that the script derives them from the data:

* n_pca = knee point of the cumulative explained-variance curve,
          capped at 50 components.
* k     = int(0.025 * n_cells) — matches the upstream tutorial.
"""
from __future__ import annotations

import argparse
import pickle
import traceback

import anndata as ad
import celloracle as co
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Sentinel value that tells the script to derive the parameter from
# the data. Kept in sync with the Nextflow module (params.n_pca_components
# and params.knn_k default to the same string).
AUTO = "auto"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-adata",   required=True)
    p.add_argument("--cluster-col",   default="leiden")
    p.add_argument("--n-pca",         default=AUTO,
                   help='Number of PCA components or "auto"')
    p.add_argument("--knn-k",         default=AUTO,
                   help='KNN k or "auto" (=> 0.025 * n_cells)')
    p.add_argument("--knn-n-jobs",    type=int,   default=4)
    p.add_argument("--grn-alpha",     type=float, default=10.0)
    p.add_argument("--filter-p",      type=float, default=0.001)
    p.add_argument("--filter-max",    type=int,   default=10000)
    p.add_argument("--output-oracle", required=True)
    p.add_argument("--output-links",  required=True)
    p.add_argument("--output-plot",   required=True)
    return p.parse_args()


def _resolve_pca_dim(n_pca_arg: str, oracle: co.Oracle) -> int:
    """Pick the number of PCA components.

    Strategy follows the upstream CellOracle tutorial: find the knee
    point of the cumulative explained-variance curve, i.e. the first
    component where the second derivative drops below ~0.002.
    """
    if str(n_pca_arg).lower() != AUTO:
        return int(n_pca_arg)

    exp_var = np.asarray(oracle.pca.explained_variance_ratio_)[:100]
    exp_var_cum = np.cumsum(exp_var)
    # second difference of the cumulative variance
    diffs = np.diff(np.diff(exp_var_cum) > 0.002)
    if diffs.size == 0 or diffs.sum() == 0:
        n_comps = int(np.argmax(exp_var_cum >= 0.90)) + 1
    else:
        n_comps = int(np.where(diffs)[0][0]) + 1
    n_comps = max(2, min(n_comps, 50))  # safety clamps
    print(f"[AUTO] PCA components selected from knee point: {n_comps}")
    return n_comps


def _resolve_knn_k(knn_k_arg: str, oracle: co.Oracle) -> int:
    """Pick the KNN k for imputation."""
    if str(knn_k_arg).lower() != AUTO:
        return int(knn_k_arg)

    n_cells = oracle.adata.shape[0]
    k = max(2, int(0.025 * n_cells))
    print(f"[AUTO] KNN k from 2.5% of {n_cells} cells: {k}")
    return k


def main() -> None:
    args = parse_args()

    # Load the pre-processed AnnData and restore the raw count layer
    # (Oracle requires un-logged counts).
    with open(args.input_adata, "rb") as fh:
        adata: ad.AnnData = pickle.load(fh)
    adata.X = adata.layers["counts"].copy()

    # Base mouse GRN distributed with CellOracle
    base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name=args.cluster_col,
        embedding_name="X_umap",
    )
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()

    # Resolve dynamic parameters
    n_comps = _resolve_pca_dim(args.n_pca, oracle)
    k       = _resolve_knn_k(args.knn_k, oracle)

    # KNN-imputation of dropouts
    oracle.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=args.knn_n_jobs,
    )

    # Cluster-wise GRN regression
    links = oracle.get_links(
        cluster_name_for_GRN_unit=args.cluster_col,
        alpha=args.grn_alpha,
        verbose_level=10,
    )
    links.filter_links(
        p=args.filter_p,
        weight="coef_abs",
        threshold_number=args.filter_max,
    )
    links.get_network_score()

    try:
        available_tfs = list(links.merged_score.index.get_level_values(0).unique())
        
        if "Mafb" in available_tfs:
            chosen_tf = "Mafb"
        elif len(available_tfs) > 0:
            chosen_tf = links.merged_score.groupby(level=0)['degree_all'].max().idxmax()
            print(f"[co_grn_inference] Using top hub TF: {chosen_tf}")
        else:
            chosen_tf = None

        if chosen_tf and hasattr(links, "plot_score_per_cluster"):
            links.plot_score_per_cluster(goi=chosen_tf, save=f"{args.output_plot}", plt_show=False)
        else:
            plt.figure()
            plt.title("Cluster network score (No TFs passed the GRN filters)")
            plt.savefig(args.output_plot, dpi=300)
            plt.close()

        cluster_ids = sorted(list(links.links_dict.keys()))
        rank_cluster = str(cluster_ids[2] if len(cluster_ids) >= 3 else cluster_ids[0])
        
        df_rank = links.merged_score[links.merged_score["cluster"] == rank_cluster]
        df_rank = df_rank.sort_values(by="degree_centrality_all", ascending=False).head(30)
        
        plt.figure(figsize=(6, 8))
        plt.scatter(df_rank["degree_centrality_all"], range(len(df_rank)))
        plt.yticks(range(len(df_rank)), df_rank.index)
        plt.xlabel("Degree Centrality (All)")
        plt.title(f"TF Rank: Cluster {rank_cluster}")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(args.output_plot.replace("grn_score.png", f"grn_rank_cluster{rank_cluster}.png"), dpi=300)
        plt.close()

        if len(cluster_ids) >= 2:
            c1, c2 = str(cluster_ids[0]), str(cluster_ids[1])
            
            df_c1 = links.merged_score[links.merged_score["cluster"] == c1]["degree_centrality_all"]
            df_c2 = links.merged_score[links.merged_score["cluster"] == c2]["degree_centrality_all"]
            
            common_tfs = df_c1.index.intersection(df_c2.index)
            x_vals = df_c1.loc[common_tfs].values
            y_vals = df_c2.loc[common_tfs].values
            
            plt.figure(figsize=(7, 6))
            plt.scatter(x_vals, y_vals, alpha=0.6, edgecolors='none', color='tab:blue')
            
            if len(common_tfs) > 0:
                combined_score = x_vals + y_vals
                threshold = np.percentile(combined_score, 99)
                for tf, x, y, score in zip(common_tfs, x_vals, y_vals, combined_score):
                    if score >= threshold:
                        plt.text(x + 0.005, y + 0.005, tf, fontsize=8)
            
            plt.xlabel(f"Degree Centrality - Cluster {c1}")
            plt.ylabel(f"Degree Centrality - Cluster {c2}")
            plt.title(f"Network Centrality Comparison: {c1} vs {c2}")
            
            max_val = max(x_vals.max() if x_vals.size > 0 else 0, y_vals.max() if y_vals.size > 0 else 0) * 1.15
            max_val = max(max_val, 0.1)
            plt.xlim(-0.01, max_val)
            plt.ylim(-0.01, max_val)
            
            plt.tight_layout()
            plt.savefig(args.output_plot.replace("grn_score.png", f"grn_comparison_{c1}vs{c2}.png"), dpi=300)
            plt.close()

    except Exception as exc:  # noqa: BLE001
        print(f"[co_grn_inference] diagnostic plot failed, generating fallback: {exc}")
        traceback.print_exc()
        plt.figure()
        plt.title("Cluster network score (Fallback due to error)")
        plt.savefig(args.output_plot, dpi=300)
        plt.close()

    # Serialise both objects — the next process needs the fitted Oracle
    # *and* the Links object to identify cluster-specific TFs.
    with open(args.output_oracle, "wb") as fh:
        pickle.dump(oracle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    with open(args.output_links, "wb") as fh:
        pickle.dump(links, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[co_grn_inference] kept {len(links.links_dict)} clusters, "
          f"{sum(len(v) for v in links.links_dict.values())} edges total")


if __name__ == "__main__":
    main()