#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

library(Seurat)
library(seuratTools)
library(infercnv)

normal_reference_path <- "~/Homo_sapiens/infercnv/reference_mat.rds"

# normal_reference_mat <- readRDS(normal_reference_path)
# 
# normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat, assay = "gene") %>%
#   RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", "."))) |> 
# 	DietSeurat(layers = "counts")

normal_seu <- readRDS(normal_reference_path)

make_seus_from_cellranger <- function(sample_path){
  # browser()
  seu_path <- path("output/seurat", paste0(path_file(sample_path), "_seu.rds"))

  print(seu_path)

  mypath <- fs::path(sample_path, "outs/filtered_feature_bc_matrix")

  count_mat <- Seurat::Read10X(mypath)

  seu <- Seurat::CreateSeuratObject(count_mat, assay = "gene") %>%
    RenameCells(new.names = str_replace(colnames(.), "-", "."))

  # saveRDS(seu, seu_path)
  print(glue::glue("saved {seu_path}"))

  return(seu)

}

cellranger_paths <-
    fs::dir_ls("output/cellranger/", glob = "*SRR*") %>%
    purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(cellranger_paths, make_seus_from_cellranger)

append_infercnv_to_seu <- function(sample_id, normal_seu, seurat_assay = "gene"){
  browser()
	
	infercnv_dir <- fs::path("output/infercnv", path_file(sample_id))
	
	infercnv_obj <- readRDS(fs::path(infercnv_dir, "/run.final.infercnv_obj"))
	
  seu_path <- fs::path("output/seurat", paste0(fs::path_file(sample_id), "_seu.rds"))

  seu <- readRDS(seu_path) |> 
  	DietSeurat(layers = "counts")

  seu_merged <- merge(seu, normal_seu, collapse = TRUE) |> 
  	# JoinLayers() |> 
  	identity()
  
  new_cell_names <- str_replace(colnames(seu_merged), "-", "\\.")
  
  seu_merged <- RenameCells(seu_merged, new.names = new_cell_names)
  
  seu_merged <- seu_merged[,colnames(seu_merged) %in% colnames(infercnv_obj@expr.data)]

  seu_merged <- infercnv::add_to_seurat(seu_merged, assay_name = seurat_assay, infercnv_dir)

  seu_w_cnv <- seu_merged[,!grepl("normal", colnames(seu_merged))]

  seu_cnv_path <- fs::path("output/seurat", paste0(path_file(sample_id), "_cnv_seu.rds"))

  seu_w_cnv0 <- seuratTools::clustering_workflow(seu_w_cnv, resolution = c(0.2, 0.4))

  # saveRDS(seu_w_cnv, seu_cnv_path)

  return(seu_w_cnv0)

}

test0 <- append_infercnv_to_seu("SRR14800534", normal_seu)



