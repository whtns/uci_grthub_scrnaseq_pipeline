#!/usr/bin/env bash
set -euo pipefail

CONTAINER="/dfs9/ucightf-lab/kstachel/singularity_containers/r_seurat_anndataR.sif"

usage() {
	cat <<'EOF'
Usage: convert_seu_to_anndata.sh --input INPUT.rds [--output OUTPUT.h5ad] [--assay ASSAY] [--x-mapping LAYER]

Convert a Seurat object saved as an RDS file to an AnnData H5AD file using
the Seurat and anndataR packages inside the Singularity container:
  /dfs9/ucightf-lab/kstachel/singularity_containers/r_seurat_anndataR.sif

Options:
  --input, -i       Path to input Seurat .rds file
  --output, -o      Path to output .h5ad file
                    Default: same basename as input with .h5ad extension
  --assay, -a       Seurat assay to export
                    Default: Seurat::DefaultAssay(seurat_obj)
  --x-mapping, -x   Seurat layer or slot to write to AnnData .X
                    Default: data
  --help, -h        Show this help text

Examples:
  ./src/convert_seu_to_anndata.sh -i output/seurat/sample.rds
  ./src/convert_seu_to_anndata.sh -i output/seurat/sample.rds -o output/scanpy/sample.h5ad
  ./src/convert_seu_to_anndata.sh -i output/seurat/sample.rds -a RNA -x counts
EOF
}

INPUT=""
OUTPUT=""
ASSAY=""
X_MAPPING="data"

while [[ $# -gt 0 ]]; do
	case "$1" in
		-i|--input)
			INPUT="$2"
			shift 2
			;;
		-o|--output)
			OUTPUT="$2"
			shift 2
			;;
		-a|--assay)
			ASSAY="$2"
			shift 2
			;;
		-x|--x-mapping)
			X_MAPPING="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Unknown argument: $1" >&2
			usage >&2
			exit 2
			;;
	esac
done

if [[ -z "$INPUT" ]]; then
	echo "Error: --input is required." >&2
	usage >&2
	exit 2
fi

if [[ ! -f "$INPUT" ]]; then
	echo "Error: input file not found: $INPUT" >&2
	exit 1
fi

if [[ -z "$OUTPUT" ]]; then
	OUTPUT="${INPUT%.*}.h5ad"
fi

mkdir -p "$(dirname "$OUTPUT")"

if ! command -v module >/dev/null 2>&1; then
	source /etc/profile >/dev/null 2>&1 || true
fi

module load singularity

if [[ ! -f "$CONTAINER" ]]; then
	echo "Error: Singularity container not found: $CONTAINER" >&2
	exit 1
fi

echo "Converting Seurat to AnnData"
echo "  Input:  $INPUT"
echo "  Output: $OUTPUT"
if [[ -n "$ASSAY" ]]; then
	echo "  Assay:  $ASSAY"
else
	echo "  Assay:  active default assay"
fi
echo "  .X:     $X_MAPPING"

SEURAT_INPUT="$(readlink -f "$INPUT")"
ANNDATA_OUTPUT="$(readlink -m "$OUTPUT")"
SEURAT_ASSAY="$ASSAY"
ANNDATA_X_MAPPING="$X_MAPPING"
export SEURAT_INPUT ANNDATA_OUTPUT SEURAT_ASSAY ANNDATA_X_MAPPING

# Build bind-mount arguments so both the input and output directories
# are accessible inside the container regardless of site auto-bind config.
INPUT_DIR="$(dirname "$SEURAT_INPUT")"
OUTPUT_DIR="$(dirname "$ANNDATA_OUTPUT")"
mkdir -p "$OUTPUT_DIR"

BIND_ARGS="${INPUT_DIR}:${INPUT_DIR}"
if [[ "$OUTPUT_DIR" != "$INPUT_DIR" ]]; then
	BIND_ARGS="${BIND_ARGS},${OUTPUT_DIR}:${OUTPUT_DIR}"
fi

singularity exec --bind "$BIND_ARGS" "$CONTAINER" Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages({
	library(Seurat)
	library(anndataR)
})

input_file <- Sys.getenv("SEURAT_INPUT")
output_file <- Sys.getenv("ANNDATA_OUTPUT")
assay_name <- Sys.getenv("SEURAT_ASSAY")
x_mapping <- Sys.getenv("ANNDATA_X_MAPPING", unset = "data")

if (!file.exists(input_file)) {
	stop("Input Seurat RDS file not found: ", input_file)
}

message("Reading Seurat object: ", input_file)
seurat_obj <- readRDS(input_file)

if (!inherits(seurat_obj, "Seurat")) {
	stop("Input object is not a Seurat object: ", input_file)
}

if (!nzchar(assay_name)) {
	assay_name <- Seurat::DefaultAssay(seurat_obj)
}

message("Using assay: ", assay_name)
message("Writing AnnData .X from: ", x_mapping)

anndataR::write_h5ad(
	seurat_obj,
	output_file,
	assay_name = assay_name,
	x_mapping = x_mapping
)

message("Wrote AnnData file: ", output_file)
EOF
