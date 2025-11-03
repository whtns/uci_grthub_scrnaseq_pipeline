#!/usr/bin/env bash
set -euo pipefail

# Generate filtering timeline plots for every sample found under output/cellranger
OUTDIR="results/qc"
mkdir -p "$OUTDIR"

for h5 in output/cellranger/*/outs/filtered_feature_bc_matrix.h5; do
	# skip if no matches
	[ -e "$h5" ] || continue
	# sample directory is two levels up from the h5 file
	sample=$(basename "$(dirname "$(dirname "$h5")")")
	out="$OUTDIR/filtering_timeline_${sample}.png"
	echo "Processing sample: $sample"
	./scripts/plot_filtering_timeline.py --input "$h5" --batch-value "$sample" --out "$out"
done

echo "All plots written to $OUTDIR"