#!/bin/bash

# Lightweight wrapper to create a CellRanger aggregation CSV from
# `output/cellranger/*/outs` and run `cellranger aggr`.
#
# This script is safe to run interactively or submit as a job. It will
# look for `molecule_info.h5` inside each sample's `outs/` directory and
# produce a CSV with the header `sample_id,molecule_h5` compatible with
# CellRanger v8+ aggregation. If a sample lacks `molecule_info.h5` it is
# skipped (but other fallbacks are attempted where sensible).

set -euo pipefail
IFS=$'\n\t'

BASE_DIR="output/cellranger"
DRY_RUN=0
ID=""
OUTDIR=""
CSVPATH=""
CELLRANGER_BIN="cellranger"

usage() {
	cat <<EOF
Usage: $(basename "$0") -i AGG_ID [-d OUTDIR] [-c CSVPATH] [-p CELLRANGER_BIN] [-n] [-h]

Create an aggregation CSV from 
  ${BASE_DIR}/*/outs
and run:
  ${CELLRANGER_BIN} aggr --id AGG_ID --csv CSVPATH

Options:
  -i AGG_ID         (required) id for the cellranger aggr run
  -d OUTDIR         directory to place CSV and to use as aggr output (default: output/cellranger_aggr_<AGG_ID>)
  -c CSVPATH        write CSV to this path (default: <OUTDIR>/aggr.csv)
  -p CELLRANGER_BIN path to cellranger executable (default: cellranger on PATH)
  -n                dry-run: only print CSV contents and the command (do not run)
  -h                show this help

Notes:
  - CSV header defaults to: sample_id,molecule_h5 (compatible with CellRanger v8).
  - The script prefers 'outs/molecule_info.h5'. If that file is missing for a
	sample it will try 'outs/filtered_feature_bc_matrix.h5' as a fallback and
	warn. If nothing suitable is found the sample is skipped.
EOF
}

while getopts ":i:d:c:p:nh" opt; do
	case "$opt" in
		i) ID="$OPTARG" ;;
		d) OUTDIR="$OPTARG" ;;
		c) CSVPATH="$OPTARG" ;;
		p) CELLRANGER_BIN="$OPTARG" ;;
		n) DRY_RUN=1 ;;
		h) usage; exit 0 ;;
		:) echo "Option -$OPTARG requires an argument." >&2; usage; exit 2 ;;
		\?) echo "Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
	esac
done

if [ -z "$ID" ]; then
	echo "Error: AGG_ID is required (-i)." >&2
	usage
	exit 2
fi

if [ -z "$OUTDIR" ]; then
	OUTDIR="output/cellranger_aggr_${ID}"
fi
mkdir -p "$OUTDIR"

if [ -z "$CSVPATH" ]; then
	CSVPATH="$OUTDIR/aggr.csv"
fi

echo "Building aggregation CSV in: $CSVPATH"

SAMPLES=()
LINES=()

if [ ! -d "$BASE_DIR" ]; then
	echo "Base directory $BASE_DIR does not exist. Nothing to aggregate." >&2
	exit 1
fi

for sample_dir in "$BASE_DIR"/*; do
	# skip non-directories
	[ -d "$sample_dir" ] || continue
	# allow cases where outputs are directly under sample_dir/outs
	outs_dir="$sample_dir/outs"
	if [ ! -d "$outs_dir" ]; then
		# sometimes cellranger outputs are nested or sample_dir may already be an outs dir
		if [ "$(basename "$sample_dir")" = "outs" ]; then
			outs_dir="$sample_dir"
			sample_root="$(dirname "$sample_dir")"
			sample_id="$(basename "$sample_root")"
		else
			echo "Warning: no outs/ in $sample_dir; skipping" >&2
			continue
		fi
	fi

	# determine library id (use sample dir name)
	if [ -z "${sample_id-}" ]; then
		sample_id="$(basename "$sample_dir")"
	fi

	# prefer molecule_info.h5
	if [ -f "$outs_dir/molecule_info.h5" ]; then
		molpath="$outs_dir/molecule_info.h5"
	elif [ -f "$outs_dir/molecule_info.h5.gz" ]; then
		molpath="$outs_dir/molecule_info.h5.gz"
	elif [ -f "$outs_dir/filtered_feature_bc_matrix.h5" ]; then
		molpath="$outs_dir/filtered_feature_bc_matrix.h5"
		echo "Note: using filtered_feature_bc_matrix.h5 for $sample_id (fallback)" >&2
	elif [ -f "$outs_dir/filtered_feature_bc_matrix.h5ad" ]; then
		molpath="$outs_dir/filtered_feature_bc_matrix.h5ad"
		echo "Note: using filtered_feature_bc_matrix.h5ad for $sample_id (fallback)" >&2
	else
		echo "Warning: no usable molecule file in $outs_dir; skipping sample $sample_id" >&2
		unset sample_id
		continue
	fi

	# write absolute path
	molpath="$(readlink -f "$molpath")"
	LINES+=("$sample_id,$molpath")
	SAMPLES+=("$sample_id")
	unset sample_id
done

if [ ${#LINES[@]} -eq 0 ]; then
	echo "No samples found to aggregate. Exiting." >&2
	exit 1
fi

# write CSV
{
	echo "sample_id,molecule_h5"
	for l in "${LINES[@]}"; do
		echo "$l"
	done
} > "$CSVPATH"

echo "Wrote CSV with ${#LINES[@]} samples:" 
for s in "${SAMPLES[@]}"; do echo " - $s"; done

CMD=("$CELLRANGER_BIN" aggr "--id" "$ID" "--csv" "$CSVPATH")

echo
echo "Prepared cellranger aggr command:"
printf ' %q' "${CMD[@]}"
echo

if [ "$DRY_RUN" -eq 1 ]; then
	echo "Dry-run: not executing cellranger aggr."; exit 0
fi

# check cellranger binary
if ! command -v "${CELLRANGER_BIN}" >/dev/null 2>&1; then
	echo "Error: ${CELLRANGER_BIN} not found on PATH. Provide -p /path/to/cellranger or install it." >&2
	exit 3
fi

echo "Running cellranger aggr (this may take a while)..."
"${CMD[@]}"

echo "cellranger aggr finished. Output is under: ${ID} (or ${OUTDIR} depending on CellRanger behavior)"

