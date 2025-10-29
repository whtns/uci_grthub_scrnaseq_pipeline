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

set -euo pipefail

# Initialize mamba/conda (if available) and activate the environment once
. ~/.mymambainit-24.3.0 2>/dev/null || true
if command -v mamba >/dev/null 2>&1; then
	mamba activate scvi-tools || true
elif command -v conda >/dev/null 2>&1; then
	conda activate scvi-tools || true
fi

SRC_NOTEBOOK="src/inspect_integrated_anndata.ipynb"

# Define destination notebooks and their corresponding combined_adata_path values
# Format: dst|path
pairs=(
	"src/inspect_integrated_anndata_adult.ipynb|../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/combined_filtered_tables_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_adult_area_1.ipynb|../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_1_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_adult_area_2.ipynb|../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_2_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_adult_area_3.ipynb|../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_3_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_adult_area_4.ipynb|../output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/area_4_harmony_integrated.h5ad"

	"src/inspect_integrated_anndata_aged.ipynb|../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/combined_filtered_tables_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_aged_area_1.ipynb|../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_1_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_aged_area_2.ipynb|../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_2_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_aged_area_3.ipynb|../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_3_harmony_integrated.h5ad"
	"src/inspect_integrated_anndata_aged_area_4.ipynb|../output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/area_4_harmony_integrated.h5ad"
)

for pair in "${pairs[@]}"; do
	dst="${pair%%|*}"
	path="${pair#*|}"
	echo "Running papermill: src=$SRC_NOTEBOOK -> dst=$dst (combined_adata_path=$path)"
	papermill "$SRC_NOTEBOOK" "$dst" -p combined_adata_path "$path"
done

echo "All papermill runs completed."