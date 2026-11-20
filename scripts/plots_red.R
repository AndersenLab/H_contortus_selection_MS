# =============================================================================
# FST selection blocks - H. contortus
# Identifies and visualizes genomic regions under selection near beta-tubulin
#
# Threshold: 99th percentile of 100 kb FST over all chromosomes (incl. X),
#            pooled over the seven Farms 08-11 comparisons.
# Blocks: runs of above-threshold windows, gaps of at most 3 windows,
#         spanning at least 10 windows.
#
# Zoomed chr I panel: 10 kb pi ratio, normalized to the genome-wide median
# (normalize_and_calculate_pi_ratio), with block boundaries from the 100 kb calls.
# =============================================================================

library(tidyverse)
library(ggtext)
library(patchwork)
library(here)
source(here("scripts", "functions.R"))

# =============================================================================
# Data Loading
# =============================================================================

fst_100kb <- load_fst_long(here("data", "fst", "100kb", "fst.csv"))

# =============================================================================
# Define Constants
# =============================================================================

target_chrom <- "hcontortus_chr1_Celeg_TT_arrow_pilon"

all_chroms <- c(
  "hcontortus_chr1_Celeg_TT_arrow_pilon",
  "hcontortus_chr2_Celeg_TT_arrow_pilon",
  "hcontortus_chr3_Celeg_TT_arrow_pilon",
  "hcontortus_chr4_Celeg_TT_arrow_pilon",
  "hcontortus_chr5_Celeg_TT_arrow_pilon",
  "hcontortus_chrX_Celeg_TT_arrow_pilon"
)

# Pairs used for threshold calculation (Farms 08-11)
threshold_pairs <- c(
  "Treatment02_Pre:Treatment08_Pre",
  "Treatment02_Pre:Treatment08_Post",
  "Treatment02_Pre:Treatment09_Pre",
  "Treatment02_Pre:Treatment09_Post",
  "Treatment02_Pre:Treatment10_Pre",
  "Treatment02_Pre:Treatment10_Post",
  "Treatment02_Pre:Treatment11_Pre"
)

# All pairs for plotting (Farms 03, 04, 08-11)
signal_pairs_all <- c(
  "Treatment02_Pre:Treatment03_Pre",
  "Treatment02_Pre:Treatment04_Pre",
  "Treatment02_Pre:Treatment08_Pre",
  "Treatment02_Pre:Treatment08_Post",
  "Treatment02_Pre:Treatment09_Pre",
  "Treatment02_Pre:Treatment09_Post",
  "Treatment02_Pre:Treatment10_Pre",
  "Treatment02_Pre:Treatment10_Post",
  "Treatment02_Pre:Treatment11_Pre"
)

# Pair label order for faceting
pair_order <- c(
  "Pre-treatment 3", "Pre-treatment 4",
  "Pre-treatment 8", "Post-treatment 8",
  "Pre-treatment 9", "Post-treatment 9",
  "Pre-treatment 10", "Post-treatment 10",
  "Pre-treatment 11"
)

MAX_GAP <- 3
MIN_SPAN <- 10

# =============================================================================
# Threshold and Block Calls
# =============================================================================

# 99th percentile threshold from Farms 08-11 (all chromosomes, incl. X)
cutoff <- quantile(
  fst_100kb %>%
    filter(pair %in% threshold_pairs, chrom %in% all_chroms) %>%
    pull(fst),
  0.99, na.rm = TRUE
)
cat("99th percentile threshold:", cutoff, "\n")

# Detect blocks for all signal pairs
all_blocks <- map_dfr(signal_pairs_all, function(p) {
  df <- fst_100kb %>% filter(chrom == target_chrom, pair == p)
  call_blocks(df, cutoff, max_gap = MAX_GAP, min_span = MIN_SPAN) %>%
    mutate(pair = p, .before = 1)
})

cat("\nDetected blocks:\n")
print(all_blocks)

# =============================================================================
# Prepare Plotting Data
# =============================================================================

plot_data <- fst_100kb %>%
  filter(pair %in% signal_pairs_all) %>%
  mutate(
    midpoint_mb = midpoint / 1e6,
    chrom_short = case_when(
      str_detect(chrom, "chr1") ~ "1",
      str_detect(chrom, "chr2") ~ "2",
      str_detect(chrom, "chr3") ~ "3",
      str_detect(chrom, "chr4") ~ "4",
      str_detect(chrom, "chr5") ~ "5",
      str_detect(chrom, "chrX") ~ "X"
    ),
    chrom_short = factor(chrom_short, levels = c("1", "2", "3", "4", "5", "X"))
  ) %>%
  filter(!is.na(chrom_short)) %>%
  left_join(all_blocks, by = "pair") %>%
  mutate(
    in_block = chrom_short == "1" & !is.na(start_mb) &
      midpoint_mb >= start_mb
