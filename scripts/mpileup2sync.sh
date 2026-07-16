#!/bin/bash
#SBATCH --job-name=mpileup2sync_merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup2sync_merge.log

set -euo pipefail

module load anaconda/2020.07
source activate /home/eschaye1/data_eande106/software/conda_envs/popoolation2

INPUT_PILEUP=${1:?"Error: Please provide input pileup file as first argument"}
OUTPUT_SYNC=${2:?"Error: Please provide output sync path as second argument"}
MIN_QUAL=${3:-20}

echo "Input pileup: $INPUT_PILEUP"
echo "Output sync: $OUTPUT_SYNC"
echo "Minimum quality: $MIN_QUAL"

java -jar /home/eschaye1/data_eande106/software/conda_envs/popoolation2/share/popoolation2-1.201-0/mpileup2sync.jar \
    --input "$INPUT_PILEUP" \
    --output "$OUTPUT_SYNC" \
    --min-qual "$MIN_QUAL" 
