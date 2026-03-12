#!/usr/bin/env Rscript
suppressPackageStartupMessages({
	library(Seurat)
	library(anndataR)
})

args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
	cat("convert_h5ad_to_seurat.R - convert an AnnData (.h5ad) file to a Seurat .rds file\n")
	cat("\nUsage:\n")
	cat("  Rscript src/convert_h5ad_to_seurat.R --data <data.h5ad> [--output <file.rds>] [--assay <assay>]\n\n")
	cat("Options:\n")
	cat("  --data  Path to input data .h5ad file (default: output/scanpy/combined_harmony_integrated.h5ad)\n")
	cat("  --output  Path to output .rds file (default: same as counts but with .rds extension)\n")
	cat("  --assay  Assay name to use after conversion (default: keep converted default assay)\n")
	invisible(NULL)
}

# defaults
data_file <- NULL
output_file <- NULL
assay_name <- NULL

if (length(args) == 0) {
	# No args provided: use existing defaults
	data_file <- "output/scanpy/combined_harmony_integrated.h5ad"
	output_file <- "output/scanpy/combined_harmony_integrated.rds"
} else {
	i <- 1
	while (i <= length(args)) {
		a <- args[i]
		if (a == "-h" || a == "--help") {
			print_usage()
			quit(status = 0)
		} else if (startsWith(a, "--data=")) {
			data_file <- sub("^--data=", "", a)
		} else if (a == "--data") {
			if ((i+1) <= length(args)) {
				data_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--data requires a value")
			}
		} else if (startsWith(a, "--output=")) {
			output_file <- sub("^--output=", "", a)
		} else if (a == "--output") {
			if ((i+1) <= length(args)) {
				output_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--output requires a value")
			}
		} else if (startsWith(a, "--assay=")) {
			assay_name <- sub("^--assay=", "", a)
		} else if (a == "--assay") {
			if ((i+1) <= length(args)) {
				assay_name <- args[i+1]
				i <- i + 1
			} else {
				stop("--assay requires a value")
			}
		} else {
			stop(paste("Unknown argument:", a))
		}
		i <- i + 1
	}

	# set defaults if missing
	if (is.null(data_file)) {
		data_file <- "output/scanpy/bp_harmony_integrated.h5ad"
	}
	if (is.null(output_file)) {
		if (grepl("\\.h5ad$", data_file, ignore.case = TRUE)) {
			output_file <- sub("\\.h5ad$", ".rds", data_file, ignore.case = TRUE)
		} else {
			output_file <- paste0(data_file, ".rds")
		}
	}
}

message("Converting: ", data_file, " -> ", output_file)

if (!file.exists(data_file)) {
	stop("Input AnnData file not found: ", data_file)
}

message("Reading AnnData with anndataR: ", data_file)
seu <- anndataR::read_h5ad(data_file, as = "Seurat")

if (!inherits(seu, "Seurat")) {
	stop("Converted object is not a Seurat object")
}

if (!is.null(assay_name) && nzchar(assay_name)) {
	if (!(assay_name %in% names(seu@assays))) {
		stop("Requested assay not found in Seurat object: ", assay_name)
	}
	Seurat::DefaultAssay(seu) <- assay_name
}

if ("batch" %in% colnames(seu[[]])) {
	Idents(seu) <- "batch"
}

saveRDS(seu, file = output_file)

message("Done. Saved Seurat object to: ", output_file)
