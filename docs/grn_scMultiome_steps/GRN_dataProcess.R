#####################################################################################
# Script Name: GRN_dataProcess.R
# Purpose: Preprocess data for the Gene Regulatory Network (GRN) pipeline.
# Inputs: Seurat object and other raw data.
# Outputs: Processed data files for downstream analysis.
#####################################################################################

# ==== Load Parameters =============================================================
# Load parameters from a YAML configuration file and set working directories.
Parameters <- "./parameters_celloracle.yaml"  # Replace with relative path
library(yaml)
var <- yaml::read_yaml(Parameters)

# Define paths based on parameters.
setwd(var$Path$Home)
path_rawdata <- paste0(var$Path$Home, var$Analysis_Name, "/raw_data/")
path_plot <- paste0(var$Path$Home, var$Analysis_Name, "/output/plots/")

# ==== Load Required Libraries =====================================================
# Load necessary libraries for data processing and analysis.
set.seed(23)
library(SeuratDisk)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(cicero)
library(EnsDb.Hsapiens.v86)

# ==== Load Data ===================================================================
# Load the Seurat object and prepare ATAC counts.
SeuObj <- readRDS(var$Path$SeuObj)
atac_counts_trans <- SeuObj@assays$ATAC$counts
rownames(atac_counts_trans) <- paste0("chr", rownames(atac_counts_trans))

Idents(SeuObj) <- "CellType"

# ==== Prepare scMultiome for scRNA ================================================
# Extract RNA data from the Seurat object and save it in H5Seurat and H5AD formats.
DefaultAssay(object = SeuObj) <- "RNA"
seu_rna <- Seurat::DietSeurat(
  SeuObj,
  counts = TRUE,  # Save raw counts to adata.layers['counts']
  data = TRUE,    # Save normalized counts to adata.X
  scale.data = FALSE,  # Avoid scaling to prevent export issues
  features = VariableFeatures(SeuObj),  # Export all genes
  assays = "RNA",
  dimreducs = c("umap_rna"),
  graphs = c("RNA_MonoOmicInteg"),
  misc = TRUE
)

# Convert factors to characters in metadata to avoid issues during export.
i <- sapply(seu_rna@meta.data, is.factor)
seu_rna@meta.data[i] <- lapply(seu_rna@meta.data[i], as.character)

# Save RNA data in H5Seurat and H5AD formats.
SaveH5Seurat(seu_rna, filename = paste0(path_rawdata, "seuobj_RNA.h5Seurat"), overwrite = TRUE)
Convert(paste0(path_rawdata, "seuobj_RNA.h5Seurat"), assay = "RNA", dest = "h5ad", overwrite = TRUE)

# ==== Prepare scMultiome for scATAC ===============================================
# Extract ATAC data from the Seurat object and save it in H5Seurat and H5AD formats.
DefaultAssay(object = SeuObj) <- "ATAC"
seu_atac <- Seurat::DietSeurat(
  SeuObj,
  counts = TRUE,  # Save raw counts to adata.layers['counts']
  data = TRUE,    # Save log1p counts to adata.X
  scale.data = FALSE,  # Avoid scaling to prevent export issues
  features = rownames(SeuObj),  # Export all peaks
  assays = "ATAC",
  dimreducs = c("umap_atac"),
  graphs = c("ATAC_MonoOmicInteg"),
  misc = TRUE
)

# Convert factors to characters in metadata to avoid issues during export.
i <- sapply(seu_atac@meta.data, is.factor)
seu_atac@meta.data[i] <- lapply(seu_atac@meta.data[i], as.character)

# Save ATAC data in H5Seurat and H5AD formats.
SaveH5Seurat(seu_atac, filename = paste0(path_rawdata, "seuobj_ATAC.h5Seurat"), overwrite = TRUE)
Convert(paste0(path_rawdata, "seuobj_ATAC.h5Seurat"), assay = "ATAC", dest = "h5ad", overwrite = TRUE)

# ==== CICERO Analysis =============================================================
# Transform the Seurat object into a CellDataSet for Cicero analysis.
DefaultAssay(object = SeuObj) <- "ATAC"
atac_cds <- as.cell_data_set(SeuObj)
atac_cds@int_colData$reducedDims@listData$UMAP <- atac_cds@int_colData$reducedDims@listData$UMAP_ATAC
atac_cds@int_colData$reducedDims@listData$UMAP_ATAC <- NULL

# Detect genes and filter out peaks with zero reads.
atac_cds <- monocle3::detect_genes(atac_cds)
atac_cds <- atac_cds[Matrix::rowSums(exprs(atac_cds)) != 0, ]

# Plot histogram of read counts for quality control.
pdf(paste0(path_plot, "Histogram_Cicero.pdf"), width = 15, height = 15)
hist(Matrix::colSums(exprs(atac_cds)))
dev.off()

# Filter cells based on read count thresholds.
max_count <- 28000
min_count <- 4000
atac_cds <- atac_cds[, Matrix::colSums(exprs(atac_cds)) >= min_count]
atac_cds <- atac_cds[, Matrix::colSums(exprs(atac_cds)) <= max_count]

# Process the Cicero CellDataSet object.
set.seed(23)
atac_cds <- detect_genes(atac_cds)
atac_cds <- estimate_size_factors(atac_cds)
atac_cds <- preprocess_cds(atac_cds, method = "LSI")

# Perform dimensional reduction using UMAP.
atac_cds <- reduce_dimension(atac_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
umap_coords <- reducedDims(atac_cds)$UMAP

# Create a Cicero object and save it.
cicero_cds <- make_cicero_cds(atac_cds, reduced_coordinates = umap_coords)
saveRDS(cicero_cds, file = paste0(path_rawdata, "cicero_cds.rds"))

# ==== Load Reference Genome Information ============================================
# Load chromosome lengths from the EnsDb.Hsapiens.v86 database.
ref <- seqlengths(EnsDb.Hsapiens.v86)
chromosome_length <- data.frame(V1 = names(ref), V2 = ref)
rownames(chromosome_length) <- 1:dim(chromosome_length)[1]
chromosome_length <- chromosome_length[nchar(chromosome_length$V1) <= 2, ]

# ==== Run Cicero ==================================================================
# Run Cicero to calculate coaccessibility scores.
conns <- run_cicero(cicero_cds, chromosome_length)

# Save Cicero results.
saveRDS(conns, file = paste0(path_rawdata, "cicero_connections.rds"))

# ==== Save Results ================================================================
# Save all peaks and Cicero connections for downstream analysis.
all_peaks <- row.names(exprs(atac_cds))
write.csv(x = all_peaks, file = paste0(path_rawdata, "all_peaks.csv"))
write.csv(x = conns, file = paste0(path_rawdata, "cicero_connections.csv"))

# Plot histogram of coaccessibility scores.
pdf(paste0(path_plot, "Histogram_Coaccess.pdf"), width = 15, height = 15)
hist(conns$coaccess[conns$coaccess != 0])
dev.off()

# Print summary statistics for coaccessibility scores.
sum(conns$coaccess == 0, na.rm = TRUE) / length(conns$coaccess)  # Fraction of zero scores.
sum(conns$coaccess > 0.1, na.rm = TRUE)  # Number of scores above 0.1.
sum(conns$coaccess == 0, na.rm = TRUE)/length(conns$coaccess) #0.225073
sum(conns$coaccess > 0.1, na.rm = TRUE) #359210

