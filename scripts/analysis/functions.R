# Shared functions for the H. contortus diversity and FST analyses.
library(tidyverse)

# --------------------------
# Define load_fst_long()
# Input: grenedalf fst.csv
# Output: long tibble with one row per window per pair (chrom, start, end, midpoint, pair, fst)

load_fst_long <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(chrom, start, end, ends_with(".fst")) %>%
    rename_with(~ str_remove(., "\\.fst$"), ends_with(".fst")) %>%
    pivot_longer(cols = -c(chrom, start, end), names_to = "pair", values_to = "fst") %>%
    mutate(midpoint = (start + end) / 2)
}
# --------------------------

# --------------------------
# Define call_blocks()
# Input: one pair's windows, FST cutoff, max gap (windows) between passing windows,
#        minimum block span (windows)
# Output: one row per block: start/end (Mb), span in windows, number of passing windows

call_blocks <- function(df, cutoff, max_gap = 3, min_span = 10) {
  df <- df %>%
    arrange(midpoint) %>%
    mutate(above = fst > cutoff)

  pass_idx <- which(df$above)
  if (length(pass_idx) == 0) {
    return(tibble(start_mb = numeric(), end_mb = numeric(),
                  span_windows = integer(), n_above = integer()))
  }

  group_id <- c(1, cumsum(diff(pass_idx) > (max_gap + 1)) + 1)

  tibble(idx = pass_idx, group_id = group_id) %>%
    group_by(group_id) %>%
    summarise(start_idx = min(idx), end_idx = max(idx), n_above = n(), .groups = "drop") %>%
    mutate(
      start_mb = df$midpoint[start_idx] / 1e6,
      end_mb = df$midpoint[end_idx] / 1e6,
      span_windows = end_idx - start_idx + 1
    ) %>%
    filter(span_windows >= min_span) %>%
    select(start_mb, end_mb, span_windows, n_above)
}
# --------------------------

# --------------------------
# Define normalize_and_calculate_pi_ratio()
# Input: Diversity values from load_diversity(), and the name of the baseline sample
# Output: Pi ratio (compared to baseline) for each non-baseline sample, raw and log10 adjusted,
#         plus the genome-wide rank of each window within that sample

normalize_and_calculate_pi_ratio <- function(data, baseline) {

  # Step 1: Normalize pi to genome-wide median (per sample)
  normalized_data <- data %>%
    group_by(sample) %>%
    mutate(pi_norm = theta_pi / median(theta_pi, na.rm = TRUE)) %>%
    ungroup()

  # Step 2: Extract baseline
  baseline_data <- normalized_data %>%
    dplyr::filter(sample == baseline) %>%
    dplyr::select(chrom, start, end, pi_norm_baseline = pi_norm)

  # Step 3: Extract non-baseline samples
  other_data <- normalized_data %>%
    dplyr::filter(sample != baseline) %>%
    dplyr::select(chrom, start, end, sample, pi_norm_sample = pi_norm)

  # Step 4: Join and calculate ratios
  ratio_data <- other_data %>%
    inner_join(baseline_data, by = c("chrom", "start", "end")) %>%
    mutate(
      pi_ratio = pi_norm_sample / pi_norm_baseline,
      pi_ratio_log10 = log10(pi_ratio)
    ) %>%
    group_by(sample) %>%
    mutate(pi_rank = rank(pi_ratio, na.last = "keep")) %>%
    ungroup() %>%
    dplyr::select(chrom, start, end, sample, pi_norm_baseline, pi_norm_sample,
                  pi_ratio, pi_ratio_log10, pi_rank)

  ratio_data
}

# --------------------------
# Define load_diversity()
# Input: Raw diversity output from grenedalf (diversity.csv)
# Output: A list of three dataframes (windows, qc, values)

load_diversity <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)

  windows <- raw %>%
    dplyr::select(chrom, start, end, starts_with("total."))

  long <- raw %>%
    dplyr::select(chrom, start, end, matches("^[^.]+\\.(missing|numeric|passed|theta_pi|theta_watterson|tajimas_d)$")) %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = c("sample", "metric"),
      names_pattern = "^(.+)\\.(missing|numeric|passed|theta_pi|theta_watterson|tajimas_d)$",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value))

  wide <- long %>%
    pivot_wider(names_from = metric, values_from = value)

  qc <- wide %>%
    dplyr::select(chrom, start, end, sample, missing, numeric, passed)

  values <- wide %>%
    dplyr::select(chrom, start, end, sample, theta_pi, theta_watterson, tajimas_d)

  list(windows = windows, qc = qc, values = values)
}
# --------------------------
