#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=merge_bams_%j.log

set -euo pipefail

# Usage: sbatch merge_bam.sh <BAM_LIST> <OUT_DIR>

module load samtools/1.15.1

BAM_LIST=${1:?"Error: Please provide BAM list as first argument"}
OUT_DIR=${2:?"Error: Please provide output directory as second argument"}

mkdir -p $OUT_DIR

# Extract unique sample prefixes
samples=$(cat $BAM_LIST | xargs -n1 basename | sed 's/_S.*//' | sort -u)

for sample in $samples; do
    echo "$(date): Merging ${sample}..."
    
    # Find all BAMs for this sample from the list
    bams=$(grep "/${sample}_S" $BAM_LIST)
    
    count=$(echo "$bams" | wc -l)
    echo "  Found ${count} BAM files:"
    echo "$bams" | sed 's/^/    /'
    
    if [[ -z "$bams" ]]; then
        echo "  WARNING: No BAMs found for ${sample}, skipping"
        continue
    fi
    
    samtools merge -@ 4 -f \
        ${OUT_DIR}/${sample}_merged.bam \
        $bams
    
    samtools index ${OUT_DIR}/${sample}_merged.bam
    
    echo "$(date): Finished ${sample}"
done

ls ${OUT_DIR}/*_merged.bam > ${OUT_DIR}/merged_bam_list.txt

echo "Done! $(wc -l < ${OUT_DIR}/merged_bam_list.txt) samples merged"
