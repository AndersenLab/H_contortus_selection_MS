library(tidyverse)
library(here)
source(here("scripts", "analysis", "functions.R"))

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

baseline  <- "Treatment02_Pre"
QUANTILE  <- 0.99   # F_ST percentile used as the selection threshold
MAX_GAP   <- 3      # windows below threshold tolerated inside a region
MIN_SPAN  <- 10     # minimum region span, in windows

chrom_order <- c(
  "hcontortus_chr1_Celeg_TT_arrow_pilon",
  "hcontortus_chr2_Celeg_TT_arrow_pilon",
  "hcontortus_chr3_Celeg_TT_arrow_pilon",
  "hcontortus_chr4_Celeg_TT_arrow_pilon",
  "hcontortus_chr5_Celeg_TT_arrow_pilon",
  "hcontortus_chrX_Celeg_TT_arrow_pilon"
)

threshold_samples <- c("Treatment08_Pre", "Treatment08_Post",
                       "Treatment09_Pre", "Treatment09_Post",
                       "Treatment10_Pre", "Treatment10_Post",
                       "Treatment11_Pre")

plot_samples <- c("Treatment03_Pre", "Treatment04_Pre",
                  "Treatment08_Pre", "Treatment08_Post",
                  "Treatment09_Pre", "Treatment09_Post",
                  "Treatment10_Pre", "Treatment10_Post",
                  "Treatment11_Pre")

threshold_pairs <- paste0(baseline, ":", threshold_samples)
plot_pairs      <- paste0(baseline, ":", plot_samples)

dir.create(here("results"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

fst_values <- load_fst_long(here("data", "fst", "100kb", "fst.csv")) %>%
  filter(chrom %in% chrom_order, !is.na(fst))

# -----------------------------------------------------------------------------
# Threshold
# -----------------------------------------------------------------------------

cutoff <- fst_values %>%
  filter(pair %in% threshold_pairs) %>%
  pull(fst) %>%
  quantile(QUANTILE, na.rm = TRUE)

writeLines(
  paste0("99th percentile F_ST threshold: ", format(cutoff, digits = 6)),
  here("results", "selected_regions_threshold.txt")
)

# -----------------------------------------------------------------------------
# Call regions and flag windows, one pair and chromosome at a time
# -----------------------------------------------------------------------------

process_pair_chrom <- function(df) {
  df <- arrange(df, midpoint)
  blocks <- call_blocks(df, cutoff, max_gap = MAX_GAP, min_span = MIN_SPAN)
  pos <- df$midpoint / 1e6
  
  in_block <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(blocks))) {
    in_block <- in_block | (pos >= blocks$start_mb[i] & pos <= blocks$end_mb[i])
  }
  
  list(
    windows = mutate(df, in_block = in_block),
    blocks  = if (nrow(blocks) > 0) {
      mutate(blocks, pair = df$pair[1], chrom = df$chrom[1])
    } else {
      NULL
    }
  )
}

results <- fst_values %>%
  filter(pair %in% plot_pairs) %>%
  group_split(pair, chrom) %>%
  map(process_pair_chrom)

# -----------------------------------------------------------------------------
# Window table
# -----------------------------------------------------------------------------

windows_out <- map_dfr(results, "windows") %>%
  mutate(
    sample      = str_extract(pair, "(?<=:).*"),
    position_mb = ((start + end) / 2) / 1e6,
    threshold   = cutoff
  ) %>%
  select(sample, pair, chrom, start, end, midpoint, position_mb,
         fst, in_block, threshold)

write_csv(windows_out, here("results", "selected_regions_windows.csv"))

# -----------------------------------------------------------------------------
# Block table
# -----------------------------------------------------------------------------

blocks_out <- map_dfr(results, "blocks") %>%
  mutate(sample = str_extract(pair, "(?<=:).*")) %>%
  select(sample, pair, chrom, start_mb, end_mb, span_windows, n_above) %>%
  arrange(sample, chrom, start_mb)

write_csv(blocks_out, here("results", "selected_regions_blocks.csv"))
