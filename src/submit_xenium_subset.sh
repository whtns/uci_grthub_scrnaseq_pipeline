#!/bin/bash
#SBATCH --job-name=xenium_combine
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
mamba activate samui

python src/subset_xenium_by_coordinates_and_combine.py -x data/output-XETG00221__0069982__Aged_cohort__20251003__181018