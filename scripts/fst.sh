#!/bin/bash
#SBATCH --job-name=grenedalf_fst
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=fst_%j.log

set -euo pipefail

SYNC=${1:?"Error: provide sync file"}
OUTDIR=${2:?"Error: provide output directory"}

module load anaconda3
source activate /home/eschaye1/data_eande106/software/conda_envs/grenedalf

mkdir -p "$OUTDIR"

echo "Sync file: $SYNC"
echo "Output directory: $OUTDIR"
echo ""

grenedalf fst \
    --sync-path "$SYNC" \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 2 \
    --filter-sample-max-read-depth 500 \
    --window-type interval \
    --window-interval-width 10000 \
    --window-interval-stride 10000 \
    --method unbiased-nei \
    --window-average-policy valid-snps \
    --pool-sizes 3000 \
    --threads 4 \
    --out-dir "$OUTDIR"

echo ""
echo "=== Done ==="
