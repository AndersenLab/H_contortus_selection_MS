#!/bin/bash
#SBATCH --job-name=grenedalf_diversity
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=diversity_%j.log

set -euo pipefail

SYNC=${1:?"Error: provide sync file"}
FORMAT=${2:?"Error: provide format (grenedalf or popoolation)"}
OUTDIR=${3:?"Error: provide output directory"}

module load anaconda3
source activate /home/eschaye1/data_eande106/software/conda_envs/grenedalf

FORMAT_FLAG=""
if [[ "$FORMAT" == "popoolation" ]]; then
    FORMAT_FLAG="--popoolation-format"
elif [[ "$FORMAT" != "grenedalf" ]]; then
    echo "Error: format must be 'grenedalf' or 'popoolation', got '$FORMAT'" >&2
    exit 1
fi
 
mkdir -p "$OUTDIR"

echo "Sync file: $SYNC"
echo "Format: $FORMAT"
echo "Output directory: $OUTDIR"
echo ""

grenedalf diversity \
    --sync-path "$SYNC" \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 2 \
    --filter-sample-max-read-depth 500 \
    --window-type interval \
    --window-interval-width 100000 \
    --window-interval-stride 100000 \
    --window-average-policy valid-loci \
    --pool-sizes 3000 \
    $FORMAT_FLAG \
    --out-dir "$OUTDIR"

echo ""
echo "=== Done ==="
