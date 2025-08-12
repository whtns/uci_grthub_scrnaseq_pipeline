#!/usr/bin/Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

library(tidyverse)
library(velocyto.R)
library(Seurat)
library(rprojroot)
library(fs)
# library(seuratTools, lib.loc = "/dataVolume/storage/rpkgs/devel_install/")
library(seuratTools)
proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

seu <- readRDS(seu_file)

# ------------------------------------------------------------------------
# IMPORTANT!!! create a loom file using the velocyto command line tools first
#line 114 cannot run in chunk; use terminal

# loom_path <- fs::path(proj_dir, "output", "velocyto", "MyTissue.loom") %>%
#   identity()

# ldat <- read.loom.matrices(loom_path)

velocity_seu <- seuratTools::velocyto_seu(seu, loom_path)

if(identical(dim(velocity_seu$gene), dim(seu$gene))){
  saveRDS(velocity_seu, str_replace(seu_file, "seurat", "velocyto"))  
}

