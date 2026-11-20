#!/bin/bash
#SBATCH --job-name=subsample_chrom
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=subsample_chrom_%A_%a.log

set -euo pipefail

module load anaconda/2020.07
source activate popoolation2

SYNC_DIR=${1:?"Error: provide sync directory"}
THRESHOLD=${2:?"Error: provide threshold"}
OUTPUT_DIR=${3:?"Error: provide output directory"}

SCRIPT=$(command -v subsample-synchronized.pl)

# Look at the sync directory to find contig files and build an array 
CHROM_FILES=($(ls "$SYNC_DIR"/*.sync | sort))

# Subtract 1 from SLURM task ID (1-based) to make it match the array (0-based) 
CHROM_FILE=${CHROM_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Extract contig name for file naming
CHROM_NAME=$(basename "$CHROM_FILE" .sync)

# Read column list from setup
COLUMNS=$(cat "$OUTPUT_DIR/columns.txt")

echo "=== Processing $CHROM_NAME ==="
echo "Input: $CHROM_FILE"
echo "Columns: $COLUMNS"
echo "Start time: $(date)"

# Take columns from samples above coverage threshold and write to temporary file
TEMP_FILE="$OUTPUT_DIR/by_contig/${CHROM_NAME}_temp.sync"
cut -f"$COLUMNS" "$CHROM_FILE" > "$TEMP_FILE"

echo "Extracted columns: $(date)"

# Subsample the sync file
OUTPUT_FILE="$OUTPUT_DIR/by_contig/${CHROM_NAME}_${THRESHOLD}x.sync"

perl "$SCRIPT" \
    --input "$TEMP_FILE" \
    --output "$OUTPUT_FILE" \
    --target-coverage "$THRESHOLD" \
    --max-coverage 2% \
    --method withoutreplace

echo "Subsampling complete: $(date)"

# Remove temporary file
rm "$TEMP_FILE"

echo "Output: $OUTPUT_FILE"
echo "End time: $(date)"
