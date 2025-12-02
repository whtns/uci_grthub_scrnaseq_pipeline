#!/usr/bin/env bash
#SBATCH --job-name=scanpy_integration
#SBATCH --output=logs/scanpy_integration_%j.out
#SBATCH --error=logs/scanpy_integration_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1


mkdir -p output/scanpy
python {params.script} \
    --filtered_matrix_dirs {input.filtered_matrix_dirs} \
    --output_prefix {params.output_prefix} \
    --min_genes {params.min_genes} \
    --min_cells {params.min_cells} \
    --n_top_genes {params.n_top_genes} \
    --batch_key {params.batch_key}