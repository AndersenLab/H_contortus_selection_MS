#!/bin/bash
#SBATCH --job-name=subsample_setup
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=subsample_setup_%j.log

set -euo pipefail

COVERAGE_FILE=${1:?"Error: provide coverage file"}
THRESHOLD=${2:?"Error: provide downsampling threshold"}
OUTPUT_DIR=${3:?"Error: provide output directory"}

echo "=== Setup ==="
echo "Coverage file: $COVERAGE_FILE"
echo "Downsampling threshold: ${THRESHOLD}x"
echo "Start time: $(date)"

mkdir -p "$OUTPUT_DIR/by_contig"

# Create sample list and column list
# Columns 1-3 are always kept (chrom, pos, ref)
# Sample columns start at 4

awk -v thresh="$THRESHOLD" -v outdir="$OUTPUT_DIR" '
BEGIN {
    FS = "\t" #Input is tab separated
    OFS = "\t" #Output is tab separated
    cols = "1,2,3" #First three columns are not sample-specific
    new_col = 4 #Column four and on are per-sample
}
NR == 1 {
    # Print header for sample list
    print "column", "ID", "timepoint", "original_coverage", "new_coverage" > outdir"/sample_list.txt"
    next
}
$2 >= thresh {
    # Value in column 2 of coverage file (specifies coverage) must be equal to or greater than the threshold provided 
    original_col = NR + 2 # Corresponding column in sync file is the row in the coverage file + 2
    cols = cols "," original_col # Include new column in the list of columns passing threshold
    print new_col, $1, $3, $2, thresh >> outdir"/sample_list.txt" # Append sample list with new column 
    new_col++ #Next sample gets next column
}
END {
    # Write column list for cut command
    print cols > outdir"/columns.txt" # Print to a file
    print "Samples passing threshold: " (new_col - 4) #Print the number of samples retained to the log 
}
' "$COVERAGE_FILE" #Input for preceeding awk command

echo ""
echo "Sample list:"
cat "$OUTPUT_DIR/sample_list.txt"
echo ""
echo "Columns to extract:"
cat "$OUTPUT_DIR/columns.txt"
echo ""
echo "End time: $(date)"
