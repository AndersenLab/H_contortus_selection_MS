#!/bin/bash
#SBATCH --job-name=mpileup2sync
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=mpileup2sync_%j.log

set -euo pipefail

INPUT_PILEUP=$1
OUTPUT_SYNC=$2
MIN_QUAL=${3:-20}

module load anaconda/2020.07
source activate popoolation2

JAR=$(find "${CONDA_PREFIX}" -name mpileup2sync.jar -print -quit)

mkdir -p "$(dirname "${OUTPUT_SYNC}")"

java -jar "${JAR}" \
    --input "${INPUT_PILEUP}" \
    --output "${OUTPUT_SYNC}" \
    --min-qual "${MIN_QUAL}"
