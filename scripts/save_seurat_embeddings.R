#!/usr/bin/Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)

embeddings_path = str_replace(seu_path, "_seu.rds", "_embeddings.csv")

seu <- readRDS(seu_path)

seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.1", "-1"))

mynb <- readRDS(nb_path)

nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

umap_embeddings <- seu@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

seurat_meta <- seu@meta.data[c("gene_snn_res.0.2", "clone_opt", "GT_opt")] %>%
  tibble::rownames_to_column("cell") %>% 
  dplyr::rename(seurat_cluster = `gene_snn_res.0.2`) %>% 
  dplyr::filter(!is.na(clone_opt))

# seurat_meta <- seu@meta.data[c("gene_snn_res.0.2")] %>%
#   tibble::rownames_to_column("cell")

umap_embeddings <- 
  umap_embeddings %>% 
  dplyr::left_join(seurat_meta, by = "cell")

write_csv(umap_embeddings, embeddings_path)

