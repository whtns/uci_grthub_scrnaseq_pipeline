#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)
library(seuratTools)
library(infercnv)

normal_reference_path <- "~/Homo_sapiens/infercnv/reference_mat.rds"

normal_reference_mat <- readRDS(normal_reference_path)
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
    RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

# make_seus_from_cellranger <- function(sample_path){
#     # browser()
#     seu_path <- path("output/seurat", paste0(path_file(sample_path), "_seu.rds"))
#
#     print(seu_path)
#
#     mypath <- fs::path(sample_path, "outs/filtered_feature_bc_matrix")
#
#     count_mat <- Seurat::Read10X(mypath)
#
#     seu <- Seurat::CreateSeuratObject(count_mat, assay = "gene") %>%
#         RenameCells(new.names = str_replace(colnames(.), "-", "."))
#
#     saveRDS(seu, seu_path)
#     print(glue::glue("saved {seu_path}"))
#
#     return(seu)
#
# }
#
# cellranger_paths <-
#     fs::dir_ls("output/cellranger/", glob = "*SRR*") %>%
#     purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))
#
# seus <- purrr::map(cellranger_paths, make_seus_from_cellranger)

clone_post_files <- dir_ls("output/numbat", regexp = "\\/SRR[0-9].*\\/clone_post_2.tsv", recurse = TRUE)

append_numbat_to_seu <- function(clone_post_file){
    # browser()
    seu_cnv_path <- path("output/seurat", paste0(path_file(path_dir(clone_post_file)), "_seu.rds"))

    seu <- readRDS(seu_cnv_path)

    clone_post <- read_tsv(clone_post_file) %>%
        dplyr::mutate(cell = str_replace_all(cell, "-", "."))  %>%
        column_to_rownames("cell")

    seu0 <- Seurat::AddMetaData(seu, clone_post)

    numbat_seu_path <- str_replace(seu_cnv_path, "_seu.rds", "_infercnv_numbat_seu.rds")

    saveRDS(seu0, numbat_seu_path)

    return(seu_cnv_path)

}

map(clone_post_files, append_numbat_to_seu)
