#!/bin/bash
#SBATCH --job-name=integrate_%j
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --error=snakemake_master.%j.err
#SBATCH --output=snakemake_master.%j.out

. ~/.mymambainit-24.3.0

mamba activate scvi-tools

papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_adult.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/combined_filtered_tables_harmony_integrated.h5ad"
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_adult_area_1.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_1_harmony_integrated.h5ad"
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_adult_area_2.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_2_harmony_integrated.h5ad"
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_adult_area_3.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_3_harmony_integrated.h5ad"
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_adult_area_4.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_4_harmony_integrated.h5ad"

papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_aged.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/combined_filtered_tables_harmony_integrated.h5ad"
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_aged_area_1.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_1_harmony_integrated.h5ad"	
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_aged_area_2.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_2_harmony_integrated.h5ad"	
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_aged_area_3.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_3_harmony_integrated.h5ad"	
papermill src/inspect_integrated_anndata.ipynb src/inspect_integrated_anndata_aged_area_4.ipynb -p combined_adata_path "../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_4_harmony_integrated.h5ad"	