#!/bin/bash
#SBATCH --job-name=indel_regions
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=indel_regions_%j.log

set -euo pipefail

module load anaconda/2020.07
source activate popoolation2

INPUT_PILEUP=${1:?"Error: Please provide input pileup file as first argument"}
OUTPUT_GTF=${2:?"Error: Please provide output GTF path as second argument"}
INDEL_WINDOW=${3:-5}
MIN_COUNT=${4:-2}

mkdir -p "$(dirname "$OUTPUT_GTF")"

echo "Input pileup: $INPUT_PILEUP"
echo "Output GTF: $OUTPUT_GTF"
echo "Indel window: $INDEL_WINDOW"
echo "Min count: $MIN_COUNT"

identify-indel-regions.pl \
    --input "$INPUT_PILEUP" \
    --output "$OUTPUT_GTF" \
    --indel-window "$INDEL_WINDOW" \
    --min-count "$MIN_COUNT"

echo "$(date): Done"
