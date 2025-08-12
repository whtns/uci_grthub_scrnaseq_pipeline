
args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

library(seuratTools)
library(sceasy)
library(reticulate)

loompy <- reticulate::import('loompy')
scanpy <- reticulate::import('scanpy')

seu <- readRDS(seu_file)

# kept_metadata <- c("nCount_RNA", "nFeature_RNA", "sample_id")

kept_metadata <- c("nCount_RNA", "nFeature_RNA", "sample_id", "batch",
									 "integrated_snn_res.0.2", "integrated_snn_res.0.4", "integrated_snn_res.0.6", "integrated_snn_res.0.8",
									 "integrated_snn_res.1.2", "integrated_snn_res.1.4", "integrated_snn_res.1.6",
									 "integrated_snn_res.1.8", "integrated_snn_res.2", "integrated_snn_res.1",
									 "S.Score", "G2M.Score", "Phase")

seu$gene@meta.data <- seu$gene@meta.data[,kept_metadata]

for(j in 1:ncol(seu$gene@meta.data)){
	if(is.factor(seu$gene@meta.data[,j]) == T){
		seu$gene@meta.data[,j] = as.character(seu$gene@meta.data[,j]) # Force the variable type to be character
		seu$gene@meta.data[,j][is.na(seu$gene@meta.data[,j])] <- "NA"
	}
	if(is.character(seu$gene@meta.data[,j]) == T){
		seu$gene@meta.data[,j][is.na(seu$gene@meta.data[,j])] <- "NA"
	}
}

anndata_file = gsub(".loom", ".h5ad", loom_file)

sceasy::convertFormat(seu$gene, from="seurat", to="anndata", outFile=anndata_file)

adata <- scanpy$read(anndata_file)
adata$write_loom(loom_file)
