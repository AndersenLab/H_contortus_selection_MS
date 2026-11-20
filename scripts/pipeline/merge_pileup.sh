#!/bin/bash
#SBATCH --job-name=merge_pileup
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=merge_pileup.log

# Usage: sbatch merge_pileup.sh <CONTIG_LIST> <PILEUP_DIR> [OUT_FILE]

set -euo pipefail

CONTIG_LIST=${1:?"Error: Please provide contig list as first argument"}
PILEUP_DIR=${2:?"Error: Please provide pileup directory as second argument"}
OUT_FILE=${3:-"${PILEUP_DIR}/merged.pileup"}

# Build array of pileup files, one per contig, in the order listed
PILEUP_FILES=()
while IFS= read -r CHROM || [ -n "${CHROM}" ]; do
    [ -z "${CHROM}" ] && continue

    f="${PILEUP_DIR}/pileup_${CHROM}.pileup"
    if [ ! -f "${f}" ]; then
        echo "Error: missing ${f}" >&2
        exit 1
    fi
    PILEUP_FILES+=("${f}")
done < "${CONTIG_LIST}"

if [ ${#PILEUP_FILES[@]} -eq 0 ]; then
    echo "Error: no contigs listed in ${CONTIG_LIST}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_FILE}")"

echo "Merging ${#PILEUP_FILES[@]} pileup files into ${OUT_FILE}"

cat "${PILEUP_FILES[@]}" > "${OUT_FILE}"

echo "Done! Lines in ${OUT_FILE}: $(wc -l < "${OUT_FILE}")"
