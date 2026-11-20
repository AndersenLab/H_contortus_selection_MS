#!/bin/bash
#SBATCH --job-name=grenedalf_fst
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=fst_%j.log

set -euo pipefail

SYNC=${1:?"Error: provide sync file"}
OUTDIR=${2:?"Error: provide output directory"}
SAMPLE_LIST=${3:?"Error: provide sample_list.txt"}
WINDOW_WIDTH=${4:-10000}
WINDOW_STRIDE=${5:-10000}
REGION=${6:-}

module load anaconda3
source activate grenedalf

SAMPLE_NAMES=$(mktemp)
trap "rm -f $SAMPLE_NAMES" EXIT
SYNC_BASE=$(basename "$SYNC" .sync)
awk -F'\t' -v base="$SYNC_BASE" 'NR>1 {print base "." NR-1 "\t" $3}' "$SAMPLE_LIST" > "$SAMPLE_NAMES"

mkdir -p "$OUTDIR"

echo "Sync file: $SYNC"
echo "Sample list: $SAMPLE_LIST"
echo "Output directory: $OUTDIR"
echo "Window width: $WINDOW_WIDTH"
echo "Window stride: $WINDOW_STRIDE"
echo "Region: ${REGION:-none (whole genome)}"

echo ""

REGION_ARGS=()
if [[ -n "$REGION" ]]; then
    REGION_ARGS=(--filter-region "$REGION")
fi

grenedalf fst \
    --sync-path "$SYNC" \
    --rename-samples-list "$SAMPLE_NAMES" \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 2 \
    --filter-sample-max-read-depth 500 \
    "${REGION_ARGS[@]}" \
    --window-type interval \
    --window-interval-width "$WINDOW_WIDTH" \
    --window-interval-stride "$WINDOW_STRIDE" \
    --method unbiased-hudson \
    --window-average-policy available-loci \
    --pool-sizes 3000 \
    --threads 4 \
    --out-dir "$OUTDIR"

echo ""
echo "=== Done ==="
