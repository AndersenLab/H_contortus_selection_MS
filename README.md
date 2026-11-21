# Subtherapeutic fenbendazole treatment rapidly selects for the F200Y variant of the isotype-1 beta-tubulin gene in Haemonchus contortus
## Abstract:
The role of variation in the isotype-1 beta-tubulin gene in benzimidazole (BZ) resistance is well established. However, resistance levels observed in the field vary among parasite populations with identical beta-tubulin variants, suggesting that additional genes likely contribute to BZ resistance. Efforts to identify uncharacterized BZ resistance genes have focused on comparisons between susceptible and resistant isolates with different geographic and demographic histories, which makes it difficult to discern whether genetic differences are associated with variation in BZ susceptibility. To identify genomic regions involved in BZ resistance, we performed a whole-genome selection scan on an experimentally selected BZ-resistant population of the strongylid parasite Haemonchus contortus. A group of 12 goats was infected with the BZ-susceptible H. contortus isolate TxPh-2011-S and treated ten times with 5% to 15% of the label dose of fenbendazole (FBZ). Pooled whole-genome sequencing (WGS) was performed before and after each treatment. We found that selection near the isotype-1 beta-tubulin locus was detectable after seven subtherapeutic treatments. Further, we found that the canonical F200Y resistance variant increased from undetectable by WGS to 92% frequency after ten subtherapeutic treatments. No additional loci displayed signatures of selection, nor were other canonical resistance variants present. The TxPh-2011-S isolate might lack additional loci that contribute to BZ resistance, or low sequencing coverage and the strong signal of isotype-1 beta-tubulin limited our power to detect small-effect loci. Our study underscores the importance of variation in the isotype-1 beta-tubulin gene in BZ resistance, and the speed with which resistance can be selected by repeated subtherapeutic treatments.

## Manuscript:
XXX

## Data available elsewhere:
* FASTQs are available at the NCBI under PRJNA1559015
* BAMs and BAM indices are available on AWS Open Data, at the path https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/H_contortus_selection_MS/
  * Per-sample BAMs and indices are named with the pattern h_contortus_<sample>_.bam and h_contortus_<sample>_.bam.bai, for example: h_contortus_sample_41
  * Key to map sample IDs to treatment timepoints is stored in the README.md

## Respository contents:

### data/
* coverage/: Per-timepoint coverage calculated by mosdepth, split by chromosome
* diversity/: Pi, theta, and Tajima's D output from grenedalf diversity
  * /100kb: Calculated in 100 kb bins across the genome
  * /10kb: Calculated in 10 kb bins across chromosome 1
* downsample/: List of samples above 10X downsampling coverage with original coverage values
* fst/: FST output from grenedalf diversity
  * /100kb: Calculated in 100 kb bins across the genome
  * /10kb: Calculated in 10 kb bins across chromosome 1
* iso-1/: Genotypes for canonical isotype-1 beta tubulin resistance codons (167, 198, and 200)
* metadata/: 
  * coverage_summary.tsv: Genome-wide mean coverage for all samples
  * timeline_clean.csv: Metadata (dates and doses) for each treatment  
