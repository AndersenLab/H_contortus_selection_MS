#!/bin/bash
#SBATCH --job-name=subsample_merge
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=subsample_merge_%j.log

set -euo pipefail

THRESHOLD=${1:?"Error: provide downsampling threshold"}
OUTPUT_DIR=${2:?"Error: provide output directory"}

echo "=== Merging ==="
echo "Downsampling threshold: ${THRESHOLD}x"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"

# Merged output file
MERGED_FILE="$OUTPUT_DIR/merged_${THRESHOLD}x.sync"

# Concatenate all subsampled contig files
cat "$OUTPUT_DIR"/by_contig/*_${THRESHOLD}x.sync > "$MERGED_FILE"

echo "Merge complete: $(date)"
echo ""
echo "Merged file: $MERGED_FILE"
echo "Merged file size: $(ls -lh "$MERGED_FILE" | awk '{print $5}')"
echo "Line count: $(wc -l < "$MERGED_FILE")"
echo ""
echo "=== Summary ==="
echo "Sample list:"
cat "$OUTPUT_DIR/sample_list.txt"
echo ""
echo "End time: $(date)"
