#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --output=coverage_%A_%a.out
#SBATCH --error=coverage_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Runs mosdepth on one BAM per array task.
#
# Usage: MOSDEPTH_BIN=/path/to/mosdepth sbatch --array=0-N coverage.sh <BAM_DIR> <OUT_DIR>

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: MOSDEPTH_BIN=/path/to/mosdepth sbatch --array=0-N coverage.sh <BAM_DIR> <OUT_DIR>"
    exit 1
fi

# Runtime parameters
BAM_DIR=$1
OUT_DIR=$2
if [ ! -d "${BAM_DIR}" ]; then
    echo "Error: BAM directory not found: ${BAM_DIR}"
    exit 1
fi

# Check mosdepth binary
if [ -z "${MOSDEPTH_BIN:-}" ]; then
    echo "Error: MOSDEPTH_BIN is not set"
    echo "Usage: MOSDEPTH_BIN=/path/to/mosdepth sbatch --array=0-N coverage.sh <BAM_DIR> <OUT_DIR>"
    exit 1
fi
if [ ! -x "${MOSDEPTH_BIN}" ]; then
    echo "Error: mosdepth not found or not executable: ${MOSDEPTH_BIN}"
    exit 1
fi
MOSDEPTH_VERSION=$("${MOSDEPTH_BIN}" --version | awk '{print $2}')
if [ "${MOSDEPTH_VERSION}" != "0.2.6" ]; then
    echo "Error: expected mosdepth 0.2.6, found ${MOSDEPTH_VERSION}"
    exit 1
fi
echo "Using mosdepth ${MOSDEPTH_VERSION} (${MOSDEPTH_BIN})"

# Create output directory
mkdir -p "${OUT_DIR}"

# Get list of BAM files
BAM_FILES=("${BAM_DIR}"/*.bam)

# Get current BAM file based on array index
BAM=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
if [ -z "${BAM}" ] || [ ! -f "${BAM}" ]; then
    echo "Error: no BAM for array task ${SLURM_ARRAY_TASK_ID} in ${BAM_DIR}"
    exit 1
fi

# Get sample name (remove path and .bam)
SAMPLE=$(basename "${BAM}" .bam)

echo "Processing: ${SAMPLE}"
echo "BAM file: ${BAM}"
echo "Output directory: ${OUT_DIR}"

# Run mosdepth
"${MOSDEPTH_BIN}" \
    --threads 4 \
    "${OUT_DIR}/${SAMPLE}" \
    "${BAM}" || { echo "Error: mosdepth failed for ${SAMPLE}"; exit 1; }

echo "Done: ${SAMPLE}"
