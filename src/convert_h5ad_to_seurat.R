library(Seurat)
library(reticulate)
library(anndata)

use_condaenv("scvi-tools", conda = "/opt/apps/mamba/24.3.0/bin/mamba")

data <- read_h5ad("output/scanpy/combined_harmony_integrated.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)
saveRDS(data,"output/seurat/combined_harmony_integrated.rds")

Seurat::Reductions(data)
