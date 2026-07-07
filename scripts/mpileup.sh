#!/bin/bash
#SBATCH --job-name=mpileup_merge
#SBATCH --array=1-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup_merge_%a.log

CONTIG_LIST=${1:?"Error: Please provide contig list as first argument"}

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_LIST")

echo "Processing chromosome: $CHROM"


samtools mpileup \
    -f /home/eschaye1/Etta/H_contortus_selection/h_contortus.PRJEB506.WBPS19.genome.fa.gz \
    -b /home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/merge/pileup/merged_bam_list.txt \
    -r "$CHROM" \
    -o /home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/merge/pileup/pileup_${CHROM}.pileup
