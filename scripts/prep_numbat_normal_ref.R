#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

# wu_normal_ref <- readRDS("output/infercnv/reference_counts.rds")

# counts <- Seurat::Read10X("/dataVolume/storage/scEiad/human_droplet_data_fetal/")
# 
# meta.data <- read_csv("/dataVolume/storage/scEiad/human_droplet_data_fetal/obs.csv") %>% 
#     column_to_rownames("index")
# 
# plae_seurat <- Seurat::CreateSeuratObject(
#     counts = counts, 
#     project = "plae_human_droplet_data",
#     assay = "RNA",
#     meta.data = meta.data
# )
# 
# saveRDS(plae_seurat, "/dataVolume/storage/scEiad/human_droplet_data_fetal_seu.rds")

collin_normal <- readRDS("~/single_cell_projects/resources/collin_et_al_proj/output/seurat/SRR13633759.rds")

prep_numbat_reference <- function(seu, meta_column){
    browser()
    normal_cell_annot <-
        seu@meta.data %>%
        tibble::rownames_to_column("cell") %>% 
    	dplyr::mutate(group = .data$CellType_predict) %>% 
        dplyr::select(cell, group) %>%
        dplyr::mutate(group = replace_na(group, "NA")) %>% 
        identity()

    ref_internal = numbat::aggregate_counts(GetAssayData(seu), normal_cell_annot)
    
}

collin_reference <- prep_numbat_reference(seu0, "CellType_predict")

plae_reference <- prep_numbat_reference(seu, "CellType")

saveRDS(plae_reference, "~/Homo_sapiens/numbat/plae_ref.rds")
