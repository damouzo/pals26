# ============================================================
# Script Name: Multiome_Integration.R
# Purpose: Integrates RNA and ATAC modalities using Seurat's multi-modal integration framework.
# Inputs: Annotated Seurat object with RNA and ATAC modalities.
# Outputs: Integrated Seurat object, UMAP visualizations, and feature plots.
# ============================================================

# ==== Argument Parsing ====
# Parse command-line arguments.
args = (commandArgs(TRUE))
if (length(args) == 0) {
    stop("No arguments supplied.")
} else {
    for (i in 1:length(args)) {
        eval(parse(text = args[[i]]))
    }
}

# ==== Settings ====
# Load parameters from the provided file.
parameters <- read.table(Parameters, sep = ":", comment.char = "#", stringsAsFactors = FALSE)
for (i in 1:nrow(parameters)) {
    varname <- trimws(parameters[i, 1])  # Remove leading/trailing spaces.
    varvalue <- trimws(parameters[i, 2])
    if (grepl(",", varvalue)) {
        varvalue <- strsplit(varvalue, ",")[[1]]
    }
    assign(varname, varvalue)
}

# Define paths for Seurat objects and plots.
path_seuobjs <- paste0(Path_Output, "archive/SeuObjs/")
dir.create(paste0(Path_Output, "plots/Multiome_Integration/"), showWarnings = FALSE)
path_plot <- paste0(Path_Output, "plots/Multiome_Integration/")

# ==== Library Loading ====
# Load required libraries.
set.seed(23)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(scCustomize))
suppressMessages(library(scNextPlot))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(BuenColors))
setwd(Path_Output)

# Define custom color palette.
my_pal <- c(jdb_palette("corona"), jdb_palette("lawhoops"), jdb_palette("brewer_spectra"),
            jdb_palette("samba_color"), jdb_palette("wolfgang_basic"))

# ==== Load Seurat Object ====
# Load the annotated Seurat object with RNA and ATAC modalities.
SeuObj <- readRDS(paste0(path_seuobjs, "_SeuObj_MonoInteg_Annot_PeakCalling.rds"))

# ==== Multi-Modal Integration ====
# Perform integration of RNA and ATAC modalities using Seurat's multi-modal framework.

# Find multi-modal neighbors using RNA and ATAC embeddings.
SeuObj <- FindMultiModalNeighbors(
    SeuObj,
    reduction.list = list("harmony_rna", "harmony_atac"),
    dims.list = list(1:ncol(Embeddings(SeuObj, "harmony_rna")), 2:ncol(Embeddings(SeuObj, "harmony_atac"))),
    modality.weight.name = c("RNA.weight", "ATAC.weight"),
    verbose = TRUE
)

# Run UMAP on the weighted nearest neighbors graph.
SeuObj <- RunUMAP(SeuObj, nn.name = "weighted.nn", assay = "RNA", reduction.name = "umap_multiome")

# Perform clustering on the weighted shared nearest neighbor (wsnn) graph.
SeuObj <- FindClusters(SeuObj, graph.name = "wsnn", resolution = 0.5)

# ==== Generate Plots ====
# Generate UMAP visualizations and feature plots for the integrated data.

# Clustering plots.
pdf(paste0(path_plot, "Multiome_Clustering_CellType.pdf"), width = 8, height = 8)
DimPlot(SeuObj, reduction = "umap_multiome", group.by = "CellType", label.size = 3, cols = my_pal)
DimPlot(SeuObj, reduction = "umap_multiome", group.by = "orig.ident", label.size = 3, cols = my_pal)
DimPlot(SeuObj, reduction = "umap_multiome", group.by = "Phase", label.size = 3, cols = my_pal)
scPlot_allres(SeuObj, clustree = TRUE, UMAPs = TRUE, clustreeGenes = NULL, reduction = "umap_multiome", Graph = "wsnn")
dev.off()

# Stemness-related feature plots.
pdf(paste0(path_plot, "Stemnes_UMAP_Multiome.pdf"), width = 8, height = 8)
FeaturePlot(SeuObj, reduction = "umap_multiome", features = "Weighted_bins", label = TRUE, label.size = 3) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

# General feature plots.
features2check <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "FRIP", "CC_Diff", "percent.mt", "percent.rb")
pdf(paste0(path_plot, "UMAP_Multiome_Features.pdf"), width = 8, height = 8)
for (feature in features2check) {
    print(FeaturePlot(SeuObj, reduction = "umap_multiome", features = feature, label = TRUE, label.size = 3) +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
}
DimPlot(SeuObj, reduction = "umap_multiome", group.by = "Phase", label = TRUE, label.size = 3, cols = my_pal)
dev.off()

# ==== Save Integrated Seurat Object ====
# Save the integrated Seurat object.
DefaultAssay(SeuObj) <- "RNA"
saveRDS(SeuObj, file = paste0(path_seuobjs, "SeuObj_MultiomeIntegr.rds"))


