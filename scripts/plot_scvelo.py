#!/usr/bin/env python

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import scvelo as scv 
import scanpy as sc
import sys
import os
import pandas as pd

h5ad_file = sys.argv[1]
plot_dir = sys.argv[2]

sample_id = os.path.basename(h5ad_file).replace("_scvelo.h5ad", "")

adata = sc.read_h5ad(h5ad_file)

clones = adata.obs.clone_opt.unique()

adata.obs_names_make_unique()

adata.obs.seurat_cluster = adata.obs.seurat_cluster.astype('category')

scv.pl.velocity_embedding_grid(adata, basis='umap', color = 'seurat_cluster', save = f'{plot_dir}/{sample_id}_embedding_grid.pdf', title = f'{sample_id}', show = False, figsize = (14,10))

# scv.pl.velocity_embedding_grid(adata, basis='umap')

# for clone in clones:
#   print(clone)
#   adata_subset = adata[adata.obs.clone_opt == clone]
#   scv.pl.velocity_embedding_grid(adata_subset, basis='umap', color = ['seurat_cluster'], save = f'{plot_dir}/{sample_id}_clone{clone}_embedding_grid.pdf', title = f'{sample_id}', show = False, figsize = (14,10))

#   # scv.pl.velocity_embedding_stream(adata, basis='umap', color = ['seurat_cluster'], save = f'{proj_dir}/results/{sample_id}_embedding_stream.pdf', title = f'{sample_id}')
#   # # scv.pl.velocity(adata, var_names=['RXRG', 'CENPF'], color = ['seurat_cluster'], save = f'{proj_dir}/results/{sample_id}_embedding_genes.pdf', title = f'{sample_id}')
# 
# 
# scv.pl.velocity_embedding_grid(adata, basis='umap', color = ['seurat_cluster'], save = f'{sample_id}_embedding_grid.pdf', title = f'{sample_id}')
# scv.pl.velocity_embedding_stream(adata, basis='umap', color = ['seurat_cluster'], save = f'{sample_id}_embedding_stream.pdf', title = f'{sample_id}')
# scv.pl.velocity(adata, var_names=['RXRG', 'CENPF'], color = ['seurat_cluster'], save = f'{sample_id}_embedding_genes.pdf', title = f'{sample_id}')
