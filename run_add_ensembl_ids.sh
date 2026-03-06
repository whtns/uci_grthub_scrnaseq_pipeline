#!/bin/bash
#
# Run add_ensembl_ids on combined_harmony_integrated.h5ad
#

set -e  # Exit on error

INPUT_FILE="output/scanpy/combined_harmony_integrated.h5ad"
SPECIES="human"

echo "Adding Ensembl gene IDs to ${INPUT_FILE}"
echo "Species: ${SPECIES}"
echo ""

# Navigate to the src directory where the script is located
cd "$(dirname "$0")/src"

python run_add_ensembl_ids.py \
    "../${INPUT_FILE}" \
    --species "${SPECIES}" \
    --var-name "gene_ids" \
    --sleep 0.05 \
    --timeout 5.0

echo ""
echo "Script completed successfully!"
