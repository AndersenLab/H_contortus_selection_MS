# Shared functions for the H. contortus diversity and FST analyses.
# Requires dplyr

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
