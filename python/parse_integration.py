# %% [markdown]
# # Integrate (batch correct) single cell datasets with scanpy and scvi-tools

# %% [markdown]
# ## load required packages

# %%
import scanpy as sc
from pathlib import Path
import scvi
import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import os


# %% [markdown]
# ## set scanpy default display settings 

# %%
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')

# %% [markdown]
# ## load parse output directories for all samples
# https://support.parsebiosciences.com/hc/en-us/articles/360052794312-Scanpy-Tutorial-65k-PBMCs
# 

# %%
# Create a Path object for the directory
directory_path = Path('/dfs9/ucightf-lab/projects/LaSpA/250716_0725Bio-04_LaSpA_Parse-WT-100K/Trailmaker_results') 

# %%
compiled_regex = re.compile(r'DGE\.mtx')
data_paths = [root for root, _, files in os.walk(directory_path) 
            if 'all-sample' not in root and any(compiled_regex.match(f) for f in files)]

# %%
data_paths

# %% [markdown]
# ## Reading in data
# 
# After reading in the data we'll perform basic filtering on our expression matrix to remove low-quality cells and uninformative genes. The parameter "min_genes" will keep cells that have at least 300 genes, and similarly, "min_cells" will keep genes that are expressed in at least 5 cells.

# %%
def read_parse(data_path, min_genes=300, min_cells=5, n_top_genes=2000, batch_key='batch'):
    """
    Performs standard Scanpy preprocessing on an AnnData object.
    
    Args:
        adata (anndata.AnnData): The AnnData object to preprocess.
        min_genes (int): Minimum number of genes expressed per cell to pass filtering.
        min_cells (int): Minimum number of cells a gene must be expressed in to pass filtering.
    """
    # The DGE_filtered folder contains the expression matrix, genes, and files
    # NOTE: split-pipe versions older than 1.1.0 used 'DGE.mtx'
    adata = sc.read_mtx(data_path + '/DGE.mtx')

    # reading in gene and cell data
    gene_data = pd.read_csv(data_path + '/all_genes.csv')
    cell_meta = pd.read_csv(data_path + '/cell_metadata.csv')
    
    # find genes with nan values and filter
    gene_data = gene_data[gene_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()
    
    # remove genes with nan values and assign gene names
    adata = adata[:,notNa]
    adata.var = gene_data
    adata.var.set_index('gene_name', inplace=True)
    adata.var.index.name = None
    # Ensure no duplicate gene symbols
    adata.var_names_make_unique()
    
    # add cell meta data to anndata object
    adata.obs = cell_meta
    adata.obs.set_index('bc_wells', inplace=True)
    adata.obs.index.name = None
    adata.obs_names_make_unique()
    
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Returns the dimensions of the expression matrix (cells, genes)
    adata.shape
    
    return adata
    

# %%
# Use glob() to get a list of all files in the directory
# The '*' wildcard matches any filename

# adatas = []
# for i, data_path in enumerate(data_paths):
#     sample_id = os.path.basename(data_path)
#     adata = read_parse(data_path)
#     adata.obs["batch"] = sample_id
#     adatas.append(adata)

# # %% [markdown]
# # ## Define a function to preprocess our anndata objects

# # %%
# def sc_preprocess(adata, n_top_genes=2000, batch_key='batch', inplace=True):
#     """
#     Performs standard Scanpy preprocessing on an AnnData object.

#     Args:
#         adata (anndata.AnnData): The AnnData object to preprocess.
#         min_genes (int): Minimum number of genes expressed per cell to pass filtering.
#         min_cells (int): Minimum number of cells a gene must be expressed in to pass filtering.
#         n_top_genes (int): Number of highly variable genes to select.
#         batch_key (str): The key in `adata.obs` that stores batch information.
#     """

#     if not inplace:
#         adata = adata.copy()
    
    
#     adata.var_names_make_unique()

#     # 2. Quality Control and Filtering
#     adata.var['mt'] = adata.var_names.str.startswith('MT-')
#     sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#     # 3. Normalization and Log-Transformation
#     adata.raw = adata
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)

#     # 4. Highly Variable Genes
#     sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)

#     sc.pp.neighbors(adata)
    
#     sc.tl.umap(adata)
    
#     # Perform clustering (if not already done)
#     sc.tl.leiden(adata, resolution=0.7)
    
#     # Find marker genes for each cluster
#     sc.tl.rank_genes_groups(adata, groupby="leiden", method='wilcoxon')
    
#     return adata

# # %% [markdown]
# # ## Run pre-processing on our anndata objects

# # %%
# # for adata in adatas:
# #     sc_preprocess(adata)

# # %% [markdown]
# # ## verify that the batch variable is appended into cell-level metadata

# # %%
# adatas[0].obs.columns

# # %% [markdown]
# # ## glue every sample together with the `concat` command
# # %%
# combined_adata = ad.concat(adatas, 
#                            join='outer',
#                            uns_merge="first")

# combined_adata.write("/dfs9/ucightf-lab/projects/LaSpA/250716_0725Bio-04_LaSpA_Parse-WT-100K/combined_anndata.h5ad")

# del(adatas)

combined_adata = sc.read_h5ad("/dfs9/ucightf-lab/projects/LaSpA/250716_0725Bio-04_LaSpA_Parse-WT-100K/combined_anndata.h5ad")   
# %%
# Ensure there are no duplicate gene symbols
combined_adata.var_names_make_unique()

# %% [markdown]
# ## 1. Quality Control and Filtering

# %%
# Assume adata is your loaded AnnData object
# Calculate QC metrics
combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

combined_adata.layers["counts"] = combined_adata.X.copy()

# Filter cells and genes based on standardized values
# ideally these values are fine-tuned by inspecting QC plots above 
sc.pp.filter_cells(combined_adata, min_genes=200)
sc.pp.filter_genes(combined_adata, min_cells=3)

# %% [markdown]
# ## 2. Normalization and Log-Transformation

# %%
# Store raw counts before normalization and transformation
combined_adata.raw = combined_adata

# Normalize to a total count of 10,000 per cell
sc.pp.normalize_total(combined_adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(combined_adata)

# %%
# Preprocessing with Scanpy
sc.pp.highly_variable_genes(
    combined_adata,
    n_top_genes=2000,
    subset=True,
    batch_key="batch",
)

# %% [markdown]
# ## establish an scvi model for integration using the 'batch' variable

# %%
scvi.model.SCVI.setup_anndata(
    combined_adata,
    batch_key="batch",
)

# %% [markdown]
# ## train the integration model 
# expect this step to take a while ~ 30 min.

# %%
model = scvi.model.SCVI(combined_adata)
model.train()

# for interactive GPU usage
# model.train(accelerator="gpu",
#             devices=-1,
#             strategy="ddp_notebook_find_unused_parameters_true")

# %% [markdown]
# ## store the integrated model within our combined object

# %%
combined_adata.obsm["X_scvi"] = model.get_latent_representation()

# %% [markdown]
# ## save this intermediate integrated data

# use scVI latent space for UMAP generation
SCVI_LATENT_KEY = "X_scvi"
sc.pp.neighbors(combined_adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(combined_adata, min_dist=0.3)

# neighbors were already computed using scVI
SCVI_CLUSTERS_KEY = "leiden_scvi"
sc.tl.leiden(combined_adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.5)

# Find marker genes for each cluster
sc.tl.rank_genes_groups(combined_adata, SCVI_CLUSTERS_KEY, method='wilcoxon')

# %%
combined_adata.write("/dfs9/ucightf-lab/projects/LaSpA/250716_0725Bio-04_LaSpA_Parse-WT-100K/combined_anndata_processed.h5ad")

