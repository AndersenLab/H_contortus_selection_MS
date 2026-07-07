#!/bin/bash
#SBATCH --job-name=mpileup_merge
#SBATCH --array=1-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup_merge_%a.log

CONTIG_LIST=${1:?"Error: Please provide contig list as first argument"}
REFERENCE=${2:?"Error: Please provide reference genome as second argument"}

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_LIST")

echo "Processing chromosome: $CHROM"
echo "Reference: $REFERENCE"


samtools mpileup \
    -f "$REFERENCE" \
    -b /home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/merge/pileup/merged_bam_list.txt \
    -r "$CHROM" \
    -o /home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/merge/pileup/pileup_${CHROM}.pileup
