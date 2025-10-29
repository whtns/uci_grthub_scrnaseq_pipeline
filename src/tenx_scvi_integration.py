#!/usr/bin/env python
# coding: utf-8
# # Integrate (batch correct) single cell datasets with scanpy and scvi-tools

# %% [markdown]
# ## load required packages

# %%
import scanpy as sc
import scanpy.external as sce
from pathlib import Path
import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import scvi
import grthub_tools as gt


# %% [markdown]
# ## set scanpy default display settings 

# %%
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')

# %% [markdown]
# ## helper functions for 10x data processing

# %%
def read_10x_data(data_path, min_genes=300, min_cells=5):
    """
    Read 10x data from filtered_feature_bc_matrix directory
    """
    # Read the 10x data
    adata = sc.read_10x_mtx(
        data_path,  # Path to the directory containing matrix.mtx file
        var_names='gene_symbols',  # use gene symbols for gene names (variables names)
        cache=True  # write a cache file for faster subsequent reading
    )
    
    # Make variable names unique (in case of duplicate gene names)
    adata.var_names_unique()
    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)  # filter out cells with too few genes
    sc.pp.filter_genes(adata, min_cells=min_cells)  # filter out genes expressed in too few cells
    
    # Add mitochondrial gene information
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    return adata

def preprocess_adata(adata, n_top_genes=2000):
    """
    Basic preprocessing for single cell data
    """
    # Save raw counts
    adata.raw = adata
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    
    # Keep only highly variable genes for downstream analysis
    adata = adata[:, adata.var.highly_variable]
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    return adata

def process_with_scvi(adata, output_prefix):
    """
    Process single cell data with scVI for batch correction and integration
    """
    # Setup scVI
    scvi.model.SCVI.setup_anndata(adata, batch_key='batch')
    
    # Train scVI model
    model = scvi.model.SCVI(adata)
    model.train()
    
    # Get the latent representation
    adata.obsm["X_scvi"] = model.get_latent_representation()
    
    # Compute neighbors and UMAP on scVI latent space
    sc.pp.neighbors(adata, use_rep="X_scvi")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    
    # Save the integrated data
    adata.write(f"{output_prefix}_integrated.h5ad")
    
    # Generate plots
    sc.pl.umap(adata, color=['batch', 'leiden'], save=f"_{output_prefix}_batch_leiden.pdf")
    
    return adata

# %% [markdown]
# ## load 10x CellRanger output directories for all samples
# Using scanpy to read 10x formatted data from CellRanger count outputs
# 

# %%

# Parse command line arguments

parser = argparse.ArgumentParser(description="Integrate single cell datasets with scanpy and scvi-tools", add_help=False)
parser.add_argument('--help', '-h', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit')
parser.add_argument('--filtered_matrix_dirs', nargs='+', help='A list cellranger output directories to process')
# parser.add_argument('--input_dir', type=str, required=True, help='Directory containing sample CellRanger output folders', default='/dfs9/ucightf-lab/projects/GreiS/250930_GreiS/output/cellranger')
parser.add_argument('--output_prefix', type=str, default="output/scanpy", help='Prefix for output files (default: input_dir name)')
parser.add_argument('--min_genes', type=int, default=300, help='Minimum genes per cell for filtering')
parser.add_argument('--min_cells', type=int, default=5, help='Minimum cells per gene for filtering')
parser.add_argument('--n_top_genes', type=int, default=2000, help='Number of highly variable genes to select')
parser.add_argument('--batch_key', type=str, default='batch', help='Batch key for integration')
args = parser.parse_args()

# directory_path = Path(args.input_dir)

# if args.output_prefix is None:
#     output_prefix = str(directory_path.name)
# else:
#     output_prefix = args.output_prefix

output_prefix = args.output_prefix

# ## Reading in data
# 
# After reading in the data we'll perform basic filtering on our expression matrix to remove low-quality cells and uninformative genes. The parameter "min_genes" will keep cells that have at least 300 genes, and similarly, "min_cells" will keep genes that are expressed in at least 5 cells.

combined_adata_path = f"{output_prefix}.h5ad"
if os.path.exists(combined_adata_path):
    print(f"{combined_adata_path} already exists. Skipping data loading and concatenation.")
    combined_adata = ad.read_h5ad(combined_adata_path)
else:
    # Look for directories containing outs/filtered_feature_bc_matrix/ with 10x files
    data_paths = []
    for filtered_matrix_dir in args.filtered_matrix_dirs:
        matrix_file = os.path.join(filtered_matrix_dir, 'matrix.mtx.gz')
        features_file = os.path.join(filtered_matrix_dir, 'features.tsv.gz')
        barcodes_file = os.path.join(filtered_matrix_dir, 'barcodes.tsv.gz')
        
        if (os.path.exists(matrix_file) and 
            os.path.exists(features_file) and 
            os.path.exists(barcodes_file)):
            data_paths.append(filtered_matrix_dir)

    print(data_paths)
    if len(data_paths) == 0:
        raise ValueError("No valid 10x filtered_feature_bc_matrix directories found.")
    adatas = []
    for i, data_path in enumerate(data_paths):
        # Extract sample ID from the path (e.g., from /path/to/01001CM_011923/outs/filtered_feature_bc_matrix)
        sample_id = os.path.basename(os.path.dirname(os.path.dirname(data_path)))
        print(f"Processing sample {sample_id} from {data_path}")
        
        # Read 10x data
        adata = read_10x_data(data_path, min_genes=args.min_genes, min_cells=args.min_cells)
        
        # Add batch information
        adata.obs[args.batch_key] = sample_id
        
        # Preprocess the data
        adata = preprocess_adata(adata, n_top_genes=args.n_top_genes)
        
        adatas.append(adata)
    
    combined_adata = ad.concat(adatas, join='outer', uns_merge="first")
    combined_adata.write(combined_adata_path)
    del adatas

print(output_prefix)

gt.integrate_w_scvi(combined_adata, output_prefix)