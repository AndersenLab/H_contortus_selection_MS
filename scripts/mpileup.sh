#!/bin/bash
#SBATCH --job-name=mpileup_merge
#SBATCH --array=1-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup_merge_%a.log

module load samtools/1.10

CONTIG_LIST=${1:?"Error: Please provide contig list as first argument"}
REFERENCE=${2:?"Error: Please provide reference genome as second argument"}
BAM_LIST=${3:?"Error: Please provide BAM list as third argument"}
OUT_DIR=${4:?"Error: Please provide output directory as fourth argument"}

mkdir -p "$OUT_DIR"

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_LIST")

echo "Processing chromosome: $CHROM"
echo "Reference: $REFERENCE"
echo "BAM list: $BAM_LIST"
echo "Output: ${OUT_DIR}/pileup_${CHROM}.pileup"


samtools mpileup \
    -f "$REFERENCE" \
    -b "$BAM_LIST" \
    -r "$CHROM" \
    -o "${OUT_DIR}/pileup_${CHROM}.pileup" \
    --ignore-RG

echo "$(date): Finished $CHROM"

