#!/bin/bash
#SBATCH --job-name=mpileup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup_%A_%a.log

# Runs samtools mpileup for one contig per array task.
#
# Usage:
#   sbatch --array=1-N mpileup.sh <CONTIG_LIST> <REFERENCE> <BAM_LIST> <OUT_DIR>
#
# Example:
#   sbatch --array=1-$(wc -l < contigs.txt) mpileup.sh contigs.txt ref.fa bam_list.txt /path/to/pileups

set -euo pipefail

if [ $# -ne 4 ]; then
    echo "Usage: sbatch --array=1-N mpileup.sh <CONTIG_LIST> <REFERENCE> <BAM_LIST> <OUT_DIR>"
    exit 1
fi

CONTIG_LIST=$1
REFERENCE=$2
BAM_LIST=$3
OUT_DIR=$4

for F in "${CONTIG_LIST}" "${REFERENCE}" "${BAM_LIST}"; do
    if [ ! -f "${F}" ]; then
        echo "Error: file not found: ${F}"
        exit 1
    fi
done

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: run this with sbatch --array, not directly"
    exit 1
fi

MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}

module load samtools/1.10

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONTIG_LIST}")
if [ -z "${CHROM}" ]; then
    echo "Error: no contig for task ${SLURM_ARRAY_TASK_ID} in ${CONTIG_LIST}"
    exit 1
fi

mkdir -p "${OUT_DIR}"
OUT="${OUT_DIR}/pileup_${CHROM}.pileup"

echo "$(date): Processing ${CHROM}"
echo "Reference: ${REFERENCE}"
echo "BAM list: ${BAM_LIST}"
echo "MAPQ >= ${MIN_MAPQ}, base quality >= ${MIN_BASEQ}"
echo "Output: ${OUT}"

if ! samtools mpileup \
        -f "${REFERENCE}" \
        -b "${BAM_LIST}" \
        -r "${CHROM}" \
        -q "${MIN_MAPQ}" \
        -Q "${MIN_BASEQ}" \
        --ignore-RG \
        -o "${OUT}"; then
    echo "Error: mpileup failed for ${CHROM}"
    rm -f "${OUT}"
    exit 1
fi

echo "$(date): Finished ${CHROM}"
