#!/bin/bash
# Combines per-sample mosdepth summaries into:
#   - one <sample>.coverage file per sample (per-contig mean depth)
#   - one summary table with genome-wide mean depth per sample
#
# Usage: bash combine_coverage.sh <MOSDEPTH_DIR> <OUT_DIR> <SUMMARY_FILE>

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: bash combine_coverage.sh <MOSDEPTH_DIR> <OUT_DIR> <SUMMARY_FILE>"
    echo "  MOSDEPTH_DIR: Directory containing .mosdepth.summary.txt files"
    echo "  OUT_DIR: Directory for individual .coverage files"
    echo "  SUMMARY_FILE: Path to output summary file"
    exit 1
fi

MOSDEPTH_DIR=$1
OUT_DIR=$2
SUMMARY_FILE=$3

if [ ! -d "${MOSDEPTH_DIR}" ]; then
    echo "Error: mosdepth directory not found: ${MOSDEPTH_DIR}"
    exit 1
fi

SUMMARIES=("${MOSDEPTH_DIR}"/*.mosdepth.summary.txt)
if [ ! -f "${SUMMARIES[0]}" ]; then
    echo "Error: no .mosdepth.summary.txt files in ${MOSDEPTH_DIR}"
    exit 1
fi

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${SUMMARY_FILE}")"

# Create summary header
printf "strain\tcoverage\n" > "${SUMMARY_FILE}"

for SUMMARY in "${SUMMARIES[@]}"; do
    SAMPLE=$(basename "${SUMMARY}" .mosdepth.summary.txt)

    # Remove the _merged suffix from BAM names so strain names are clean
    STRAIN=${SAMPLE%_merged}

    echo "Processing: ${STRAIN}"

    # Individual per-contig coverage file
    awk -v STRAIN="${STRAIN}" 'BEGIN{N=0} \
        { if (NR > 1){N=N+1; CHROMS[N]=$1; COVERAGE[$1]=$4;}} \
        END{printf "strain\tcoverage"; for (I=1;I<N;I++) printf "\t%s", CHROMS[I]; printf "\n"; \
            printf "%s\t%f", STRAIN, COVERAGE[CHROMS[N]]; for (I=1;I<N;I++) printf "\t%s", COVERAGE[CHROMS[I]]; \
            printf "\n";}' "${SUMMARY}" > "${OUT_DIR}/${STRAIN}.coverage"

    # Add genome-wide mean (the "total" row) to the summary
    tail -1 "${OUT_DIR}/${STRAIN}.coverage" | cut -f1,2 >> "${SUMMARY_FILE}"
done

echo "Done!"
echo "Individual files: ${OUT_DIR}/"
echo "Summary file: ${SUMMARY_FILE}"
