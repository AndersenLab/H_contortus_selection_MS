# Pipeline scripts

### merge_bam.sh
Merges per-lane BAMs into one BAM per sample, then indexes merged BAMs

### coverage.sh
Calculates coverage for each sample (BAM)

### combine_coverage.sh
Combines per-sample mosdepth summaries (`*.mosdepth.summary.txt`) into per-sample per-contig mean depth files and one genome-wide coverage table

### mpileup.sh
Runs samtools mpileup on BAMs to make a pileup for each contig

### merge_pileup.sh
Concatenates per-contig pileups into one pileup file, in the order given by a contig list

### mark_indels.sh
Identifies genomic regions near indels from a pileup file

### mpileup2sync.sh
Converts a pileup to a sync file

### filter_indels.sh
Removes regions in close proximity to indels

### split_sync.sh
Splits a sync file into one file per contig for parallelization

### subsample_setup.sh
Selects samples whose genome-wide mean coverage meets a threshold to be used by the subsampling steps

### subsample_chrom.sh
Subsamples each contig's sync file to a target coverage

### subsample_merge.sh
Concatenates per-contig subsampled sync files into one genome-wide subsampled sync file and prints the sample list

### subsample_wrapper.sh
Submits the subsampling pipeline as three sequential jobs: setup, subsampling, and merge

### diversity.sh
Runs grenedalf diversity in windows on a subsampled sync file

### fst.sh
Runs grenedalf fst in windows on a subsampled sync file
