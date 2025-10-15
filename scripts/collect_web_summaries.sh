#!/usr/bin/env bash
set -euo pipefail

# collect_web_summaries.sh
# Find all web_summary.html files produced by cellranger and copy them
# into a single directory for easy viewing.
#
# Usage: collect_web_summaries.sh [CELLRANGER_ROOT] [DEST_DIR]
# Defaults: CELLRANGER_ROOT=output/cellranger
#           DEST_DIR=output/web_summaries

CELLR_ROOT="${1:-output/cellranger}"
DEST_DIR="${2:-output/web_summaries}"

usage(){
    cat <<EOF
Usage: $0 [CELLRANGER_ROOT] [DEST_DIR]

Finds files matching:
  CELLRANGER_ROOT/<sample>/outs/web_summary.html
and copies them to DEST_DIR as <sample>_web_summary.html

Defaults:
  CELLRANGER_ROOT=${CELLR_ROOT}
  DEST_DIR=${DEST_DIR}
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

mkdir -p "$DEST_DIR"

found=0
# Use find to locate web_summary.html files. Print0 to be safe with spaces.
while IFS= read -r -d '' src; do
    found=$((found+1))
    # Expect path like: .../cellranger/<sample>/outs/web_summary.html
    # sample is the parent directory of 'outs'
    sample_dir=$(dirname "$(dirname "${src}")")
    sample=$(basename "${sample_dir}")
    dest="$DEST_DIR/${sample}_web_summary.html"
    cp -p "${src}" "${dest}"
    echo "Copied: ${src} -> ${dest}"
done < <(find "${CELLR_ROOT}" -maxdepth 3 -type f -name "web_summary.html" -print0)

if [[ ${found} -eq 0 ]]; then
    echo "No web_summary.html files found under ${CELLR_ROOT}" >&2
    exit 1
fi

# Mark completion
touch "${DEST_DIR}/.done"
echo "Collected ${found} web_summary.html file(s) into ${DEST_DIR}"

exit 0
