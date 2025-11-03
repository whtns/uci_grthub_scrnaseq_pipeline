#!/usr/bin/env python

import scanpy as sc
import sys
import os
import pandas as pd
import scvelo as scv

anndata_file = sys.argv[1]
loom_file = sys.argv[2]
seurat_embedding_path = sys.argv[3]


# anndata_file = "output/scanpy/SRR14800540.h5ad"
# loom_file = "output/velocyto/SRR14800540.loom"
# seurat_embedding_path = "output/seurat/SRR14800540_embeddings.csv"

sample_id = os.path.basename(loom_file).replace(".loom", "")

processed_anndata_file = anndata_file.replace(".h5ad", "_scvelo.h5ad")

seurat_embedding = pd.read_csv(seurat_embedding_path)
seurat_embedding.cell = seurat_embedding.cell.str.replace(".1", "-1")

adata = scv.read(anndata_file, cache=True)

seurat_embedding = seurat_embedding[seurat_embedding.cell.isin(adata.obs.index)]
adata = adata[adata.obs.index.isin(seurat_embedding.cell)].copy()

ldata = scv.read(loom_file, cache=True)

scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

adata = scv.utils.merge(adata, ldata)

# scv.pp.filter_and_normalize(adata)
# 
# scv.utils.cleanup(adata)
# 
# scv.pp.moments(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.moments(adata, method='hnsw')
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.umap(adata)

seurat_clusters = seurat_embedding[["seurat_cluster"]]
seurat_clusters.index = adata.obs.index
adata.obs["seurat_cluster"] = seurat_clusters

myclones = seurat_embedding[["clone_opt"]]
myclones.index = adata.obs.index
adata.obs["clone_opt"] = myclones

seurat_embedding = seurat_embedding[["UMAP_1", "UMAP_2"]].to_numpy()
adata.obsm["X_umap"] = seurat_embedding

scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

adata.write_h5ad(processed_anndata_file)

