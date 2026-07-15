#!/bin/bash
#SBATCH --job-name=merge_pileup
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=merge_pileup.log

source activate /home/eschaye1/data_eande106/software/conda_envs/popoolation2

CONTIG_LIST=${1:?"Error: Please provide contig list as first argument"}
PILEUP_DIR=${2:?"Error: Please provide pileup directory as second argument"}
OUT_FILE=${3:-"${PILEUP_DIR}/final_merged.pileup"}

set -euo pipefail

# Build array of pileup files, one per contig, in the order listed
PILEUP_FILES=()
while IFS= read -r CHROM; do
    f="${PILEUP_DIR}/pileup_${CHROM}.pileup"
    if [ -f "$f" ]; then
        PILEUP_FILES+=("$f")
    else
        echo "WARNING: missing $f" >&2
    fi
done < "$CONTIG_LIST"

echo "Merging ${#PILEUP_FILES[@]} pileup files into $OUT_FILE"

cat "${PILEUP_FILES[@]}" > "$OUT_FILE"

echo "Done! Lines in ${OUT_FILE}: $(wc -l < "$OUT_FILE")"
