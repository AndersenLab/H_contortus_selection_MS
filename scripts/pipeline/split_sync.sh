#!/bin/bash
#SBATCH --job-name=split_sync
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=split_sync_%j.log

set -euo pipefail

INPUT_SYNC=${1:?"Error: Please provide input sync file as first argument"}
OUTPUT_DIR=${2:?"Error: Please provide output directory  as second argument"}

mkdir -p "$OUTPUT_DIR"

echo "Input file : $INPUT_SYNC" 
echo "Start time: $(date)"

# Print contig names
echo "Chromosomes in file:"
cut -f1 "$INPUT_SYNC" | sort -u

# Split by chromosome
awk '{print > "'$OUTPUT_DIR'/"$1".sync"}' "$INPUT_SYNC"

echo "End time: $(date)"
echo "Split complete. Files created:"
ls -lh "$OUTPUT_DIR"/

