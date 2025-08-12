args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tximport)
  library('tidyverse')
  library('fs')
  library('rprojroot')
  library(seuratTools)
})

# outrds = "output/seurat/salmon_seu.rds"
organism = "Homo_sapiens"
proj_dir = "."

# print(stringtiedir)
print(proj_dir)
print(outrds)
# print(legacy_seu_rds)
print(organism)
print(getwd())

organism = c("Mus_musculus" = "mouse", "Homo_sapiens" = "human")[organism]

proj_dir = rprojroot::find_root(criterion = rprojroot::has_file_pattern("*.Rproj"))

# debug(load_counts_by_tximport)
txi_features <- load_counts_by_tximport(proj_dir, type = "stringtie")

tpm_meta <- seuratTools::load_meta(proj_dir)

seu <- seu_from_tximport(txi_features, tpm_meta)

seu <- seuratTools::clustering_workflow(seu, organism = organism)
saveRDS(seu, file = unfiltered_seu_rds)

# legacy_seus <- seuratTools::clustering_workflow(feature_seus, organism = organism, legacy_settings = TRUE)
# saveRDS(legacy_seus, file = legacy_seu_rds)

sessionInfo()
date()

