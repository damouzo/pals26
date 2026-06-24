#####################################################################################
# Script Name: GRN_analysis.py
# Purpose: Perform Python-based analysis for the GRN pipeline.
# Inputs: Processed data files and R analysis results.
# Outputs: GRN results and additional data.
#####################################################################################

# ==== Import Libraries ============================================================
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shutil
import pickle
import scanpy as sc
import seaborn as sns
import os, sys, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
from matplotlib.backends.backend_pdf import PdfPages
from gimmemotifs.motif import read_motifs
import subprocess

# Set global plotting parameters
plt.rcParams['figure.figsize'] = (15, 7)
plt.rcParams["savefig.dpi"] = 600

# ==== Load Parameters ============================================================
# Load analysis parameters from a YAML file
parameters_file = "./parameters_celloracle.yaml"  # Replace with relative path
with open(parameters_file, 'r') as file:
    var = yaml.safe_load(file)

# Define paths based on parameters
os.chdir(var['Path']['Home'])
path_rawdata = os.path.join(var['Path']['Home'], 'raw_data')
path_output = os.path.join(var['Path']['Home'], 'output')
path_plot = os.path.join(var['Path']['Home'], 'output', 'plots')
path_genome_ref = os.path.join(var['Path']['genome_ref'])

# ==== Data Preprocessing =========================================================
# Load scATAC-seq data
adata = sc.read_h5ad(os.path.join(path_rawdata, 'object_ATAC.h5ad'))
adata.layers['counts'] = adata.raw.X.copy()
adata.layers['norm'] = adata.X.copy()
peaks = adata.var.index.to_numpy()

# Load Cicero coaccessibility scores
cicero_connections = pd.read_csv(os.path.join(path_rawdata, "cicero_connections.csv"), index_col=0)

# Format peaks and Cicero connections
peaks = np.array(['chr' + peak.replace('-', '_') for peak in peaks], dtype=object)

def modify_columns(value):
    value = str(value).replace("-", "_")
    if not value.startswith("chr"):
        value = "chr" + value
    return value

cicero_connections[['Peak1', 'Peak2']] = cicero_connections[['Peak1', 'Peak2']].applymap(modify_columns)

# ==== Annotate Transcription Start Sites (TSSs) ==================================
# Annotate TSSs using CellOracle's motif analysis module
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")

# Integrate TSS info with Cicero connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)

# ==== Filter Peaks ===============================================================
# Filter peaks based on coaccessibility scores
peak = integrated[integrated.coaccess >= 1]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

# Save filtered peaks to a CSV file
peak.to_csv(os.path.join(path_output, "processed_peak_file.csv"))

# ==== Motif Analysis =============================================================
# Load processed peaks
peaks = pd.read_csv(os.path.join(path_output, "processed_peak_file.csv"), index_col=0)
ref_genome = "hg38"
peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=path_genome_ref)

# Create TFinfo object and scan for TF binding motifs
tfi = ma.TFinfo(peak_data_frame=peaks, ref_genome=ref_genome, genomes_dir=path_genome_ref)

# Load custom motif database
from gimmemotifs.motif import MotifConfig
config = MotifConfig()
motif_dir = config.get_motif_dir()
motif_db_option = "gimme.vertebrate.v5.0"
motif_db_file = motif_db_option + ".pfm"
path = os.path.join(motif_dir, motif_db_file)
motifs = read_motifs(path)

# Perform motif scanning
tfi.scan(motifs=motifs, verbose=True)

# Save TFinfo object
tfi.to_hdf5(file_path=os.path.join(path_output, "test1.celloracle.tfinfo"))

# ==== Base GRN Construction ======================================================
# Filter motifs and create base GRN
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

# Save base GRN as a DataFrame
df = tfi.to_dataframe()
df.to_parquet(os.path.join(path_output, f"base_GRN_dataframe_{motif_db_option}.parquet"))

# ==== Oracle Object Creation =====================================================
# Load scRNA-seq data
adata = sc.read_h5ad(os.path.join(path_rawdata, 'object_RNA.h5ad'))
adata.layers['counts'] = adata.raw.X.copy()
adata.layers['norm'] = adata.X.copy()

# Instantiate Oracle object and load data
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name=var['Info4GRN']['clustering'],
    embedding_name=var['Info4GRN']['umap']
)

# Load base GRN into Oracle object
base_GRN = pd.read_parquet(os.path.join(path_output, f"base_GRN_dataframe_{motif_db_option}.parquet"))
oracle.import_TF_data(TF_info_matrix=base_GRN)

# ==== GRN Calculation ============================================================
# Perform PCA and KNN imputation
oracle.perform_PCA()
n_comps = min(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_) > 0.002)))[0][0], 50)
k = int(0.025 * oracle.adata.shape[0])
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)

# Save Oracle object
oracle.to_hdf5(os.path.join(path_output, f"OracleHSC_data_{motif_db_option}.celloracle.oracle"))

# Calculate GRNs
links = oracle.get_links(cluster_name_for_GRN_unit=var['Info4GRN']['clustering'], alpha=10, verbose_level=10)

# ==== Export GRNs ================================================================
# Save GRNs as CSV and Excel files
cluster = var['Info4GRN']['clust2study']
links.links_dict[cluster].to_csv(os.path.join(path_output, f"raw_{motif_db_option}_GRN_for_{cluster}_.csv"))

output_file = os.path.join(path_output, f"raw_{motif_db_option}_GRN_clusters.xlsx")
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    for cluster, df in links.links_dict.items():
        df.to_excel(writer, sheet_name=cluster, index=False)

# Save Links object
links.to_hdf5(file_path=os.path.join(path_output, f"links_{motif_db_option}.celloracle.links"))

# ==== Network Analysis ===========================================================
# Filter network edges and calculate network scores
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
links.get_network_score()
links.merged_score.to_excel(os.path.join(path_output, "links_merged_score.xlsx"), index=True)

# Save filtered Links object
links.to_hdf5(file_path=os.path.join(path_output, f"selected_links_{motif_db_option}.celloracle.links"))

# ==== Visualization ==============================================================
# Plot network scores and comparisons
links.plot_scores_as_rank(cluster=links.cluster[0], n_gene=30, save=os.path.join(path_plot, f"Network_ranked_score"))
links.plot_score_comparison_2D(
    value="degree_centrality_all",
    cluster1="CelltypeA", cluster2="CelltypeB",
    percentile=98, save=os.path.join(path_plot, f"score_comparison")
)
plt.clf()

# -------------------------- degree centrality in 
links.plot_score_comparison_2D(
    value="degree_centrality_in",
    cluster1="CelltypeA", cluster2="CelltypeB",
    percentile=98, save=os.path.join(path_plot, f"score_comparison")
)
plt.clf()
    
# -------------------------- degree centrality out
links.plot_score_comparison_2D(
    value="degree_centrality_out",
    cluster1="CelltypeA", cluster2="CelltypeB",
    percentile=98, save=os.path.join(path_plot, f"score_comparison")
)

  # KNN imputation
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)

  # Save and Load ------
oracle.to_hdf5(os.path.join(path_output,"OracleHSC_data_{motif_db_option}.celloracle.oracle"))



# GRN calculation --------------------------------------------------
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
# Visualice the Clustering in the dimensional Reduction 
adata.obsm["X_umap"] = adata.obsm[umap].copy()
plot1 = sc.pl.umap(adata, color=clustering)
pp = PdfPages(os.path.join(path_plot,"Check_Clusterin_DimRed.pdf"))
pp.savefig(plot1)
pp.close()

# Get GRNs --------
links = oracle.get_links(cluster_name_for_GRN_unit=clustering, alpha=10, verbose_level=10) # This step may take some time.(~30 minutes)


# Export GRNs --------
cluster= var['Info4GRN']['clust2study'] 
links.links_dict.keys()
links.links_dict[cluster]
  # Save as csv
links.links_dict[cluster].to_csv(os.path.join(path_output, f"raw_{motif_db_option}_GRN_for_{cluster}_.csv")) 


# Usar ExcelWriter para escribir m�ltiples DataFrames en hojas separadas
output_file = os.path.join(path_output, f"raw_{motif_db_option}_GRN_clusters.xlsx")
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    for cluster, df in links.links_dict.items():
        df.to_excel(writer, sheet_name=cluster, index=False)



# Change order --------
  # Show the contents of pallete
links.palette

  # Save Links object.
links.to_hdf5(file_path= os.path.join(path_output, f"links_{motif_db_option}.celloracle.links"))



# Network preprocessing ----------------------------------------------------
  ## Filter network edges
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

  ## Degree distribution 
plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True, save=os.path.join(path_plot, f"Distributions"))
plt.rcParams["figure.figsize"] = [6, 4.5]

  ## Calculate netowrk score 
links.get_network_score()
links.merged_score.head()
links.merged_score.to_excel(os.path.join(path_output, "links_merged_score.xlsx"), index=True)

  ## Save Links object.
links.to_hdf5(file_path=os.path.join(path_output,f"selected_links_{motif_db_option}.celloracle.links"))


# Network analysis; Network score for each gene -----------------------------
  # Network score in each cluster, Check cluster name
links.cluster
  # Visualize top n-th genes that have high scores.
links.plot_scores_as_rank(cluster=links.cluster[0], n_gene=30, save=os.path.join(path_plot, f"Network_ranked_score"))

# Network score comparison between two clusters
plt.figure(figsize=(4,6)) 
plt.rcParams.update({'font.size': 10}) 
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
links.plot_score_comparison_2D(value="degree_centrality_all",
                               cluster1="celltypeA", cluster2="CelltypeB", 
                               percentile=98, save=os.path.join(path_plot, f"score_comparison"))
plt.clf()
                               
