#!/bin/bash
#SBATCH --job-name=merge_ab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=merge_ab_%j.log

module load samtools/1.15.1

MERGED_DIR=${1:?"Error: Please provide merged BAM directory as first argument"}

cd $MERGED_DIR

# Merge 47a + 47b
echo "$(date): Merging 47a + 47b..."
samtools merge -@ 4 47_merged.bam 47a_merged.bam 47b_merged.bam
samtools index 47_merged.bam
echo "$(date): Finished 47"

# Merge Ext_1504a + Ext_1504b
echo "$(date): Merging Ext_1504a + Ext_1504b..."
samtools merge -@ 4 Ext_1504_merged.bam Ext_1504a_merged.bam Ext_1504b_merged.bam
samtools index Ext_1504_merged.bam
echo "$(date): Finished Ext_1504"

# Remove intermediate files
echo "$(date): Removing intermediate files..."
rm 47a_merged.bam* 47b_merged.bam* Ext_1504a_merged.bam* Ext_1504b_merged.bam*

# Update the merged_bam_list.txt
ls *_merged.bam > merged_bam_list.txt

echo "$(date): Done! Final sample count: $(wc -l < merged_bam_list.txt)"
