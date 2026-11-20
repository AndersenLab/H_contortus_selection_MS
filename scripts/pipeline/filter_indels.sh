#!/bin/bash
#SBATCH --job-name=filter_indels
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=filter_indels_%j.log

set -euo pipefail

module load anaconda/2020.07
source activate popoolation2

INPUT_SYNC=${1:?"Error: Please provide input sync file as first argument"}
INPUT_GTF=${2:?"Error: Please provide input GTF file as second argument"}
OUTPUT_SYNC=${3:?"Error: Please provide output sync path as third argument"}

mkdir -p "$(dirname "$OUTPUT_SYNC")"

echo "Input sync: $INPUT_SYNC"
echo "Input GTF: $INPUT_GTF"
echo "Output sync: $OUTPUT_SYNC"

filter-sync-by-gtf.pl \
    --input "$INPUT_SYNC" \
    --gtf "$INPUT_GTF" \
    --output "$OUTPUT_SYNC"

echo "$(date): Done"


