#!/usr/bin Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}


## ----load-packages------------------------------------------------------------
library(glue)
library(fs)
library(rprojroot)
library(tidyverse)
library(seuratTools)
library(infercnv)

## -----------------------------------------------------------------------------

matrix_dir = fs::path_dir(matrix_file)

count_mat <- Seurat::Read10X(matrix_dir) %>% 
	identity()

# ref_counts <- readRDS("../output/infercnv/reference_counts.rds")

ref_counts <- readRDS(normal_reference_mat) %>% 
	identity()

combine_matrices <- function(exp_mat, reference_mat){
	exp_names <- rownames(exp_mat)
	reference_names <- rownames(reference_mat)
	shared_names <- intersect(exp_names, reference_names)
	
	exp_mat <- exp_mat[shared_names,]
	reference_mat <- reference_mat[shared_names,]
	
	colnames(reference_mat) <- paste0("normal_", colnames(reference_mat))
	
	common_mat <- cbind(exp_mat, reference_mat)
	
	colnames(common_mat) <- str_replace(colnames(common_mat), "-", ".")
	
	return(common_mat)
	
}

combined_count_mat <- combine_matrices(count_mat, ref_counts)

## -----------------------------------------------------------------------------

# create annotations files 
annotations <- 
  colnames(combined_count_mat) %>%
	as_tibble_col(column_name = "sample_id") %>%
	dplyr::mutate(annotation = ifelse(str_detect(sample_id, "normal"), "normal", "rb_tumor")) %>%
	identity()

# annotations_path <- fs::path(proj_dir, "output/infercnv/annotations_file.tsv") 

write.table(annotations, annotations_path, sep = "\t", row.names = FALSE, col.names = FALSE)

gene_order_file="~/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.gene_pos.txt" %>% 
  read_tsv(col_names = c("symbol", "seqnames", "start", "end")) %>% 
  dplyr::filter(seqnames %in% c(1:22, "X", "Y", "MT")) %>%
  identity()

combined_count_mat <- combined_count_mat[rownames(combined_count_mat) %in% gene_order_file$symbol,]

# counts_path <- fs::path(proj_dir, "output/infercnv/gene_matrix.tsv")


## -----------------------------------------------------------------------------

infercnv_obj  <- CreateInfercnvObject(raw_counts_matrix=combined_count_mat, annotations_file=annotations_path, delim = "\t", gene_order_file="~/Homo_sapiens/grch38_tran/Homo_sapiens.GRCh38.87.gene_pos.txt", ref_group_names = "normal")
	

infercnv_obj  <- infercnv::run(infercnv_obj,
                               analysis_mode='subclusters',
                               cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir = out_dir,
                               cluster_by_groups= TRUE,
                               denoise = TRUE,
                               HMM = TRUE,
                               debug = TRUE,
															 num_threads = as.integer(threads))

# infercnv_obj <- readRDS(paste0(out_dir, "/run.final.infercnv_obj"))
# 
# # -----------------------------------------------------------------------------
# plot_cnv(infercnv_obj,
#          # plot_chr_scale = TRUE,
# 				 out_dir = out_dir,
#          output_format = "png",
#          write_expr_matrix = FALSE)

