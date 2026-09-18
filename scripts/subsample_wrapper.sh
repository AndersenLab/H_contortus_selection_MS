#!/bin/bash

set -euo pipefail

SYNC_DIR=${1:?"Error: provide sync directory"}
COVERAGE_FILE=${2:?"Error: provide coverage file"}
THRESHOLD=${3:?"Error: provide threshold"}
OUTPUT_DIR=${4:?"Error: provide output directory"}

# Find directory where wrapper script lives to infer location of other scripts
SCRIPT_DIR="/home/eschaye1/Etta/H_contortus_selection_paper/scripts"

# Count contigs for array job
NUM_CONTIGS=$(ls "$SYNC_DIR"/*.sync | wc -l)

echo "Sync directory: $SYNC_DIR"
echo "Coverage file: $COVERAGE_FILE"
echo "Downsampling threshold: ${THRESHOLD}x"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Submit setup job
SETUP_JOB=$(sbatch --parsable \
    "$SCRIPT_DIR/subsample_setup.sh" \
    "$COVERAGE_FILE" \
    "$THRESHOLD" \
    "$OUTPUT_DIR")
echo "Submitted setup job: $SETUP_JOB"

# Submit downsampling job (after setup is complete)
ARRAY_JOB=$(sbatch --parsable \
    --dependency=afterok:$SETUP_JOB \
    --array=1-${NUM_CONTIGS} \
    "$SCRIPT_DIR/subsample_chrom.sh" \
    "$SYNC_DIR" \
    "$THRESHOLD" \
    "$OUTPUT_DIR")
echo "Submitted array job: $ARRAY_JOB"

# Submit merge job (after downsampling is complete)
MERGE_JOB=$(sbatch --parsable \
    --dependency=afterok:$ARRAY_JOB \
    "$SCRIPT_DIR/subsample_merge.sh" \
    "$THRESHOLD" \
    "$OUTPUT_DIR")
echo "Submitted merge job: $MERGE_JOB"

echo ""
echo "=== All jobs submitted ==="
