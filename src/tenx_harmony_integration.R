#!/usr/bin/env Rscript

# TenX/CellRanger Harmony Integration with Seurat
# Integrates 10x/CellRanger single-cell data using Harmony

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Integrate tenx data using Harmony and Seurat")
parser$add_argument("--input_dir", type = "character", default = "output/cellranger", required = FALSE, help = "Directory containing 10x/CellRanger sample output folders")
parser$add_argument("--output_prefix", type = "character", default = "cellranger_comb", help = "Prefix for output files")
parser$add_argument("--min_genes", type = "integer", default = 300, help = "Minimum genes per cell")
parser$add_argument("--min_cells", type = "integer", default = 5, help = "Minimum cells per gene")
parser$add_argument("--n_top_genes", type = "integer", default = 2000, help = "Number of variable features")
parser$add_argument("--batch_key", type = "character", default = "batch", help = "Batch variable name")
parser$add_argument("--output_dir", type = "character", default = "output/seurat", required = FALSE, help = "Output directory")
parser$add_argument("--ncores", type = "integer", default = 4, help = "Number of cores")
parser$add_argument("--harmony_theta", type = "double", default = 2, help = "Harmony theta parameter")
parser$add_argument("--harmony_dims", type = "integer", default = 30, help = "Number of Harmony dimensions")
parser$add_argument("--cluster_resolution", type = "double", default = 0.5, help = "Clustering resolution")
parser$add_argument("--globals_max_size", type = "integer", default = 24e3, help = "Clustering resolution")
parser$add_argument("--max_genes", type = "integer", default = 8000, help = "Maximum genes per cell (filter)")
parser$add_argument("--max_percent_mt", type = "double", default = 5, help = "Maximum mitochondrial percentage per cell (filter)")
parser$add_argument("--groupvar", type = "character", default = "type", help = "Column name in metadata to use for grouping (group.by)")
parser$add_argument("--group1", type = "character", default = "Medino_rooster_3NP5", help = "Value in groupvar for ident.1")
parser$add_argument("--group2", type = "character", default = "TC_med_4_7P5", help = "Value in groupvar for ident.2")

args <- parser$parse_args()
groupvar <- args$groupvar
group1 <- args$group1
group2 <- args$group2

# Check required arguments
if (is.null(args$input_dir)) {
  stop("input_dir is required")
}
if (is.null(args$output_dir)) {
  stop("output_dir is required")
}

# Set up parallel processing
library(future)
plan("multicore", workers = args$ncores)
options(future.globals.maxSize = args$globals_max_size * 1024^2)  # 24GB

cat("Starting TenX/CellRanger Harmony integration...\n")
cat("Input directory:", args$input_dir, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("Using", args$ncores, "cores\n")

# Function to read 10x / CellRanger output
read_tenx_data <- function(data_dir, min_cells = 5, min_genes = 200) {
  cat("Reading 10x/CellRanger data from:", data_dir, "\n")

  # Expect a sample folder that contains an `outs/` directory
  outs_dir <- file.path(data_dir, "outs")
  if (!dir.exists(outs_dir)) {
    stop("outs directory not found in: ", data_dir)
  }

  # Common CellRanger locations: either a folder `outs/filtered_feature_bc_matrix/`
  # or an HDF5 file `outs/filtered_feature_bc_matrix.h5`.
  ff_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  h5_file <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")

  counts <- NULL
  if (file.exists(h5_file)) {
    cat("Found HDF5 matrix:", h5_file, "\n")
    counts <- Read10X_h5(h5_file)
  } else if (dir.exists(ff_dir)) {
    cat("Found matrix folder:", ff_dir, "\n")
    counts <- Read10X(ff_dir)
  } else {
      stop("No filtered_feature_bc_matrix folder or HDF5 found in: ", outs_dir)
  }

  # Read10X may return a list (e.g., when multiple assays are present)
  if (is.list(counts) && length(counts) == 1) {
    counts <- counts[[1]]
  }

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = min_cells,
    min.features = min_genes,
    project = basename(data_dir)
  )

  return(seurat_obj)
}

# Find all sample directories with Parse output
sample_dirs <- list.dirs(args$input_dir, recursive = FALSE)
sample_dirs <- sample_dirs[sapply(sample_dirs, function(x) {
  # detect CellRanger/10x outputs: outs/filtered_feature_bc_matrix folder or HDF5
  dir.exists(file.path(x, "outs", "filtered_feature_bc_matrix")) ||
    file.exists(file.path(x, "outs", "filtered_feature_bc_matrix.h5"))
})]

# sample_dirs <- sapply(sample_dirs, function(x) {
#   file.path(x, "DGE_filtered")
# })

if (length(sample_dirs) == 0) {
  stop("No Parse output directories found in: ", args$input_dir)
}

cat("Found", length(sample_dirs), "sample directories\n")

# Read all samples
seurat_list <- list()
for (i in seq_along(sample_dirs)) {
  sample_name <- basename(sample_dirs[i])
  cat("Processing sample:", sample_name, "\n")

  seurat_obj <- read_tenx_data(sample_dirs[i], min_genes =  args$min_genes, min_cells = args$min_cells)
  seurat_obj[[]][[args$batch_key]] <- sample_name
  seurat_list[[sample_name]] <- seurat_obj
}

# Merge all samples if multiple
if (length(seurat_list) > 1) {
  cat("Merging", length(seurat_list), "samples...\n")
  combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], 
                          add.cell.ids = names(seurat_list))
} else {
  combined_seurat <- seurat_list[[1]]
}

cat("Combined object dimensions:", dim(combined_seurat), "\n")

# Add mitochondrial gene percentage
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")

# Filter cells by maximum genes and mitochondrial percent
cat("Filtering cells by nFeature_RNA <=", args$max_genes, "and percent.mt <", args$max_percent_mt, "%\n")
before_cells <- ncol(combined_seurat)
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA <= args$max_genes & percent.mt < args$max_percent_mt)
after_cells <- ncol(combined_seurat)
cat("Filtered cells: kept", after_cells, "of", before_cells, "cells\n")

# Normalize and find variable features
cat("Normalizing data...\n")
combined_seurat <- NormalizeData(combined_seurat)

cat("Finding variable features...\n")
combined_seurat <- FindVariableFeatures(combined_seurat, 
                                       selection.method = "vst",
                                       nfeatures = args$n_top_genes)

# Scale data
cat("Scaling data...\n")
combined_seurat <- ScaleData(combined_seurat, 
                            features = VariableFeatures(combined_seurat))

# Run PCA
cat("Running PCA...\n")
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(combined_seurat))

# Run Harmony integration
cat("Running Harmony integration...\n")
combined_seurat <- RunHarmony(combined_seurat, 
                             group.by.vars = args$batch_key,
                             reduction.use = "pca",
                             dims.use = 1:args$harmony_dims,
                             theta = args$harmony_theta,
                             plot_convergence = FALSE)

# Run UMAP on Harmony embeddings
cat("Running UMAP...\n")
combined_seurat <- RunUMAP(combined_seurat, 
                          reduction = "harmony", 
                          dims = 1:args$harmony_dims)

# Find clusters
cat("Finding clusters...\n")
combined_seurat <- FindNeighbors(combined_seurat, 
                                reduction = "harmony", 
                                dims = 1:args$harmony_dims)
combined_seurat <- FindClusters(combined_seurat, resolution = args$cluster_resolution)

combined_seurat <- JoinLayers(combined_seurat)

Idents(combined_seurat) <- combined_seurat$seurat_clusters

combined_seurat@misc$markers <- FindAllMarkers(combined_seurat, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

clusters <- unique(Idents(combined_seurat))
combined_seurat@misc$diffex <- lapply(clusters, function(cl) {
  sub <- subset(combined_seurat, idents = cl)
  FindMarkers(sub,
              ident.1 = group1,
              ident.2 = group2,
              group.by = groupvar,
              min.pct = 0.25,
              logfc.threshold = 0.25)
})
names(combined_seurat@misc$diffex) <- clusters

# Save integrated Seurat object
output_file <- file.path(args$output_dir, paste0(args$output_prefix, "_harmony_integrated.rds"))
cat("Saving integrated object to:", output_file, "\n")
saveRDS(combined_seurat, file = output_file)

# Save embeddings
embeddings_file <- file.path(args$output_dir, paste0(args$output_prefix, "_harmony_embeddings.csv"))
cat("Saving embeddings to:", embeddings_file, "\n")

embeddings_df <- data.frame(
  cell_id = rownames(combined_seurat@meta.data),
  combined_seurat@meta.data,
  combined_seurat@reductions$umap@cell.embeddings,
  combined_seurat@reductions$harmony@cell.embeddings[, 1:10]  # First 10 Harmony dimensions
)
write.csv(embeddings_df, file = embeddings_file, row.names = FALSE)

# Create plots
cat("Creating plots...\n")
plots_file <- file.path(args$output_dir, paste0(args$output_prefix, "_harmony_plots.pdf"))

pdf(plots_file, width = 12, height = 8)

# UMAP by batch
p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = args$batch_key) +
  ggtitle("UMAP colored by batch")
print(p1)

# UMAP by clusters
p2 <- DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP colored by clusters")
print(p2)

# Feature plots for QC metrics
p3 <- FeaturePlot(combined_seurat, features = "nFeature_RNA") +
  ggtitle("Number of features per cell")
print(p3)

p4 <- FeaturePlot(combined_seurat, features = "nCount_RNA") +
  ggtitle("Total UMI counts per cell")
print(p4)

p5 <- FeaturePlot(combined_seurat, features = "percent.mt") +
  ggtitle("Mitochondrial gene percentage")
print(p5)

# Harmony convergence plot if available
if ("harmony" %in% names(combined_seurat@tools)) {
  if (!is.null(combined_seurat@tools$harmony$plot_convergence)) {
    print(combined_seurat@tools$harmony$plot_convergence)
  }
}

dev.off()

cat("Harmony integration completed successfully!\n")
cat("Output files:\n")
cat("  - Seurat object:", output_file, "\n")
cat("  - Embeddings:", embeddings_file, "\n")
cat("  - Plots:", plots_file, "\n")