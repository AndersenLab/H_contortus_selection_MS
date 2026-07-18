#!/bin/bash
#SBATCH --job-name=grenedalf_diversity
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=diversity_%j.log

set -euo pipefail

SYNC=${1:?"Error: provide sync file"}
OUTDIR=${2:?"Error: provide output directory"}
SAMPLE_LIST=${3:?"Error: provide sample_list.txt"}

module load anaconda3
source activate /home/eschaye1/data_eande106/software/conda_envs/grenedalf


SAMPLE_NAMES=$(mktemp)
trap "rm -f $SAMPLE_NAMES" EXIT
SYNC_BASE=$(basename "$SYNC" .sync)
awk -F'\t' -v base="$SYNC_BASE" 'NR>1 {print base "." NR-1 "\t" $3}' "$SAMPLE_LIST" > "$SAMPLE_NAMES"

mkdir -p "$OUTDIR"

echo "Sync file: $SYNC"
echo "Sample list: $SAMPLE_LIST"
echo "Output directory: $OUTDIR"
echo ""

grenedalf diversity \
    --sync-path "$SYNC" \
    --rename-samples-list "$SAMPLE_NAMES" \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 2 \
    --filter-sample-max-read-depth 500 \
    --window-type interval \
    --window-interval-width 100000 \
    --window-interval-stride 100000 \
    --window-average-policy valid-loci \
    --pool-sizes 3000 \
    --out-dir "$OUTDIR"

echo ""
echo "=== Done ==="
