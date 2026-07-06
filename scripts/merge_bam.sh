#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=merge_bams_%j.log


module load samtools/1.15.1

BAM_LIST="/home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/pileup/bam_list.txt"
OUT_DIR="/home/eschaye1/Etta/H_contortus_selection/popoolation2_analyses/merge"

mkdir -p $OUT_DIR

# Extract unique sample prefixes
samples=$(cat $BAM_LIST | xargs -n1 basename | sed 's/_L00[12]\.bam//' | sort -u)

for sample in $samples; do
    echo "$(date): Merging ${sample}..."
    
    L001=$(grep "${sample}_L001" $BAM_LIST)
    L002=$(grep "${sample}_L002" $BAM_LIST)
    
    if [[ -z "$L001" || -z "$L002" ]]; then
        echo "WARNING: Could not find both lanes for ${sample}"
        continue
    fi
    
    samtools merge -@ 4 -f \
        ${OUT_DIR}/${sample}_merged.bam \
        $L001 $L002
    
    samtools index ${OUT_DIR}/${sample}_merged.bam
    
    echo "$(date): Finished ${sample}"
done

ls ${OUT_DIR}/*_merged.bam > ${OUT_DIR}/merged_bam_list.txt

echo "Done! $(wc -l < ${OUT_DIR}/merged_bam_list.txt) samples merged"
