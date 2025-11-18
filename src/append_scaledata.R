#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript src/append_scaledata.R PATH_TO_SCALED_DATA_DIRECTORY")
}

data_dir <- args[1]
if (!dir.exists(data_dir)) {
    stop(paste("Directory not found:", data_dir))
}

# Load the scaled data matrix
scaled_data <- readMM(file.path(data_dir, "scaled_data.mtx"))

# Load cell and gene names
cell_barcodes <- read.csv(file.path(data_dir, "scaled_data_barcodes.csv"))$Barcode
gene_names <- read.csv(file.path(data_dir, "scaled_data_genes.csv"))$GeneName

# Set row and column names for the matrix
rownames(scaled_data) <- gene_names
colnames(scaled_data) <- cell_barcodes

# Create a Seurat object using the scaled data
# Note: If you have raw counts or normalized data, you'd typically create the Seurat object with that
# and then use SetAssayData for the scaled data.
# For directly creating with scaled data, ensure you understand the implications for downstream analyses
# that might expect raw counts or normalized data in the 'counts' or 'data' slot.

seurat_obj <- readRDS("output/Seurat5Shiny/seurat5.rds")

# Determine target assay (use default assay if present)
assay_name <- tryCatch(DefaultAssay(seurat_obj), error = function(e) NULL)
if (is.null(assay_name) || assay_name == "") assay_name <- NULL

# Try to set the scaled data using the newer 'layer' argument first.
# Use `try(..., silent=TRUE)` (no `finally` arg) and only convert sparse->dense
# if it is required to succeed, to avoid large memory allocations when possible.
res <- try(SetAssayData(object = seurat_obj, assay = assay_name, layer = 'scale.data', new.data = scaled_data), silent = TRUE)
if (!inherits(res, "try-error")) {
    seurat_obj <- res
} else {
    message("layer-based SetAssayData failed, falling back to slot-based attempt.")
    # Attempt slot-based call without converting if scaled_data is sparse
    if (inherits(scaled_data, "sparseMatrix") || inherits(scaled_data, "dgTMatrix") || inherits(scaled_data, "dgCMatrix")) {
        rownames(scaled_data) <- rownames(GetAssayData(seurat_obj, layer = "data"))
        colnames(scaled_data) <- colnames(GetAssayData(seurat_obj, layer = "data"))
        res2 <- try(SetAssayData(object = seurat_obj, assay = assay_name, slot = 'scale.data', new.data = scaled_data), silent = TRUE)
        if (!inherits(res2, "try-error")) {
            seurat_obj <- res2
        } else {
            message("Setting sparse 'scale.data' failed: attempting dense conversion. This may require a lot of memory.")
            seurat_obj <- SetAssayData(object = seurat_obj, layer = "scale.data", as.matrix(scaled_data), assay = assay_name)
        }
    } else {
        # scaled_data already dense; set directly
        seurat_obj <- SetAssayData(object = seurat_obj, assay = assay_name, slot = 'scale.data', new.data = scaled_data)
    }
}

# Save the modified Seurat object (overwrite original)
seurat_path <- file.path("output/Seurat5Shiny/seurat5_mod.rds")
saveRDS(seurat_obj, seurat_path)
cat("Wrote updated Seurat object to:", seurat_path, "\n")