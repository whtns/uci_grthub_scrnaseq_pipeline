#!/bin/bash
#SBATCH --job-name=integrate_%j
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --error=logs/harmony_integration%j.err
#SBATCH --output=logs/harmony_integration%j.out

# Fail fast: exit on error, treat unset vars as errors, and make pipelines fail on any step
set -euo pipefail

# Resolve the directory this script lives in so relative paths inside the repo work
# even when the script is invoked from another working directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default runtime parameters (can be overridden via CLI flags below)
# - MIN_GENES: minimum number of genes per cell used by the integration script
# - RUN_INTEGRATION: toggle to skip running the Python integration step
# - INSPECT_NOTEBOOK: notebook used by papermill to inspect results
MIN_GENES=200
RUN_INTEGRATION=true
INSPECT_NOTEBOOK="src/inspect_integrated_anndata.ipynb"

# Argument parsing: supports the following options and then one or more
# `output_prefix` positional arguments. Examples of flags:
#   --min_genes N        Set minimum genes filter passed to the Python script
#   --no-integration      Skip execution of the Python integration step
#   --notebook PATH       Use a different papermill notebook for inspection
# After flags, provide one or more `output_prefix` values (paths where results live).
# Example invocation:
#   ./src/submit_harmony_integration.sh --min_genes 500 output/prefix1 output/prefix2
while [[ $# -gt 0 ]]; do
	case "$1" in
		--min_genes|--min-genes)
			if [[ -z "${2:-}" ]]; then
				echo "Error: value required for $1" >&2
				exit 2
			fi
			MIN_GENES="$2"
			shift 2
			;;
        --no-integration)
            # Skip the Python integration step; still run downstream inspection/convert
            RUN_INTEGRATION=false
            shift
            ;;
		--notebook|--inspect-notebook|--inspect_notebook)
			if [[ -z "${2:-}" ]]; then
				echo "Error: value required for $1" >&2
				exit 2
			fi
			INSPECT_NOTEBOOK="$2"
			shift 2
			;;
		--help|-h)
			# Print short usage and exit
			echo "Usage: $0 [--min_genes N] output_prefix [output_prefix ...]"
			echo
			echo "Runs src/harmony_integration.py for each provided output_prefix."
			echo "Optional flags: --inspect-notebook PATH (default: src/inspect_integrated_anndata.ipynb)"
			echo
			echo "Example prefixes:" 
			echo "  output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/control"
			exit 0
			;;
		-* )
			echo "Unknown option: $1" >&2
			exit 2
			;;
		* )
			break
			;;
	esac
done

if [[ $# -lt 1 ]]; then
	echo "Error: at least one output_prefix is required." >&2
	echo "Usage: $0 [--min_genes N] output_prefix [output_prefix ...]" >&2
	exit 2
fi

. ~/.mymambainit-24.3.0
mamba activate scvi-tools

# Loop over each provided `output_prefix` and run the integration + inspection steps.
# The Python integration step can be skipped with `--no-integration` so that
# only the report-generation (papermill/nbconvert) runs.
for output_prefix in "$@"; do
		if [[ "$RUN_INTEGRATION" == true ]]; then
			echo "---- Running integration for: $output_prefix (min_genes=$MIN_GENES) ----"
			python src/harmony_integration.py \
				--output_prefix "$output_prefix" \
				--min_genes "$MIN_GENES"
		else
			echo "Skipping python step for: $output_prefix (--no-integration was set)"
		fi
    
    slug=$(basename "$output_prefix")
    DIRNAME=$(dirname "$output_prefix")
	
	# Run papermill and nbconvert but don't let failures stop the entire job.
	# Execute papermill to render an inspection notebook for the integrated AnnData
	# (writes to "$DIRNAME/inspect_integrated_anndata_${slug}.ipynb"). Failures are
	# non-fatal: they are caught and logged so processing continues for other prefixes.
	if papermill "$INSPECT_NOTEBOOK" \
		"$DIRNAME/inspect_integrated_anndata_${slug}.ipynb" \
		-p combined_adata_path "$DIRNAME/${slug}_harmony_integrated.h5ad"; then
		echo "Papermill succeeded for ${slug}"
	else
		echo "Warning: papermill failed for ${slug}" >&2
	fi

	if jupyter nbconvert --to webpdf --allow-chromium-download --no-input \
		"$DIRNAME/inspect_integrated_anndata_${slug}.ipynb"; then
		echo "nbconvert succeeded for ${slug}"
	else
		echo "Warning: nbconvert failed for ${slug}" >&2
	fi

	echo "---- Finished: $output_prefix ----"
done
