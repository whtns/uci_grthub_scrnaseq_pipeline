#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
usage <- function() {
	cat("Usage: Rscript src/normalize_and_scale_seurat.R PATH_TO_SEURAT_OBJECT\n")
	cat("Accepts saved Seurat objects via readRDS() or .RData/.rda files containing a Seurat object.\n")
}

if (length(args) < 1) {
	usage()
	stop("Missing path to Seurat object.")
}

path <- args[1]
if (!file.exists(path)) stop(paste("File not found:", path))

seurat_obj <- NULL
# Try readRDS first (common for single-object saves)
try({
	obj <- readRDS(path)
	if (!is.null(obj) && "Seurat" %in% class(obj)) {
		seurat_obj <- obj
	}
}, silent = TRUE)

if (is.null(seurat_obj)) {
	# Try loading into a new environment (for .RData/.rda files)
	env <- new.env()
	nms <- try(load(path, envir = env), silent = TRUE)
	if (!inherits(nms, "try-error")) {
		for (nm in nms) {
			val <- env[[nm]]
			if (!is.null(val) && "Seurat" %in% class(val)) {
				seurat_obj <- val
				break
			}
		}
	}
}

if (is.null(seurat_obj)) stop("No Seurat object found in file.")

cat(sprintf("Loaded Seurat object: %s cells x %s features\n", ncol(seurat_obj), nrow(seurat_obj)))

# Normalize, find variable features, and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Save processed object next to original file with suffix
out_dir <- dirname(path)
base <- tools::file_path_sans_ext(basename(path))
out_path <- file.path(out_dir, paste0(base, "_normalized_scaled.rds"))
saveRDS(seurat_obj, out_path)
cat("Saved normalized and scaled Seurat object to:", out_path, "\n")
