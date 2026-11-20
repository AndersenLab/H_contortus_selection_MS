# =============================================================================
# H. contortus FST Selection Analysis - Final Script with Parallelization
# =============================================================================

library(tidyverse)
library(furrr)

# Set up parallel processing
plan(multisession, workers = 4)  # Adjust based on your CPU cores

# =============================================================================
# Load data
# =============================================================================

load_fst <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(chrom, start, end, ends_with(".fst")) %>%
    rename_with(~ str_remove(., "\\.fst$"), ends_with(".fst")) %>%
    pivot_longer(cols = -c(chrom, start, end), names_to = "pair", values_to = "fst") %>%
    mutate(midpoint = (start + end) / 2)
}

fst_100kb <- load_fst("~/Etta/H_contortus_selection/fst/genome_100kb/fst.csv")

# =============================================================================
# Constants
# =============================================================================

target_chrom <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
beta_tubulin_pos <- 7.025e6

autosomes_I_V <- c(
  "hcontortus_chr1_Celeg_TT_arrow_pilon",
  "hcontortus_chr2_Celeg_TT_arrow_pilon",
  "hcontortus_chr3_Celeg_TT_arrow_pilon",
  "hcontortus_chr4_Celeg_TT_arrow_pilon",
  "hcontortus_chr5_Celeg_TT_arrow_pilon"
)

signal_pairs <- c(
  "Treatment02_Pre:Treatment08_Pre",
  "Treatment02_Pre:Treatment08_Post",
  "Treatment02_Pre:Treatment09_Pre",
  "Treatment02_Pre:Treatment09_Post",
  "Treatment02_Pre:Treatment10_Pre",
  "Treatment02_Pre:Treatment10_Post",
  "Treatment02_Pre:Treatment11_Pre"
)

# =============================================================================
# Calculate threshold - 99th percentile, no exclusion
# =============================================================================

bg_data <- fst_100kb %>%
  filter(pair %in% signal_pairs, chrom %in% autosomes_I_V)

cutoff <- quantile(bg_data$fst, 0.99, na.rm = TRUE)
cat("99th percentile threshold:", round(cutoff, 4), "\n")

# =============================================================================
# Block-finding function
# =============================================================================

find_all_blocks <- function(df, cutoff, max_gap = 3) {
  df <- df %>% arrange(midpoint) %>% mutate(pass = fst > cutoff)
  pass_idx <- which(df$pass)
  
  if (length(pass_idx) == 0) {
    return(tibble(start_mb = numeric(), end_mb = numeric(), 
                  span_windows = integer(), n_pass = integer()))
  }
  
  gaps <- diff(pass_idx)
  group_id <- c(1, cumsum(gaps > (max_gap + 1)) + 1)
  
  tibble(idx = pass_idx, group_id = group_id) %>%
    group_by(group_id) %>%
    summarise(start_idx = min(idx), end_idx = max(idx), n_pass = n(), .groups = "drop") %>%
    mutate(
      start_mb = df$midpoint[start_idx] / 1e6,
      end_mb = df$midpoint[end_idx] / 1e6,
      span_windows = end_idx - start_idx + 1
    ) %>%
    select(start_mb, end_mb, span_windows, n_pass)
}

# =============================================================================
# Detect blocks
# =============================================================================

results <- map_dfr(signal_pairs, function(p) {
  blocks <- fst_100kb %>%
    filter(pair == p, chrom == target_chrom) %>%
    find_all_blocks(cutoff, max_gap = 3) %>%
    filter(span_windows >= 10)
  
  if (nrow(blocks) > 0) {
    blocks %>% 
      mutate(pair = p, width_mb = round(end_mb - start_mb, 2)) %>%
      select(pair, start_mb, end_mb, width_mb, span_windows, n_pass)
  } else {
    tibble(pair = p, start_mb = NA, end_mb = NA, width_mb = NA, 
           span_windows = NA, n_pass = NA)
  }
})

cat("\n=== Block Detection Results ===\n")
print(results)

# =============================================================================
# Permutation test (parallelized) - FIXED: genome-wide shuffle
# =============================================================================

# Pre-filter data for speed
fst_filtered <- fst_100kb %>%
  filter(pair %in% signal_pairs, chrom %in% autosomes_I_V) %>%
  arrange(pair, chrom, midpoint)

# Observed max span
obs_max <- max(results$span_windows, na.rm = TRUE)
cat("\nObserved max span:", obs_max, "windows\n")

# Run permutation test
n_perm <- 2000
cat("Running permutation test (n =", n_perm, ") with", availableCores(), "cores...\n")

set.seed(42)
start_time <- Sys.time()

null_max_spans <- future_map_dbl(1:n_perm, function(i) {
  
  # For each permutation: shuffle FST values ACROSS the entire genome
  # within each pair (breaks spatial structure entirely)
  
  null_blocks <- map_dfr(signal_pairs, function(p) {
    
    pair_data <- fst_filtered %>% filter(pair == p)
    
    # Shuffle all FST values for this pair (across all chromosomes)
    pair_data$fst <- sample(pair_data$fst)
    
    # Now look for blocks on each chromosome
    map_dfr(autosomes_I_V, function(chr) {
      chr_data <- pair_data %>% filter(chrom == chr)
      
      if (nrow(chr_data) == 0) return(tibble())
      
      find_all_blocks(chr_data, cutoff, max_gap = 3) %>%
        filter(span_windows >= 10)
    })
  })
  
  if (nrow(null_blocks) > 0) max(null_blocks$span_windows) else 0
  
}, .options = furrr_options(seed = TRUE), .progress = TRUE)

end_time <- Sys.time()
cat("Completed in", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# Calculate p-value
p_value <- mean(null_max_spans >= obs_max)

cat("\n=== Permutation Test Results ===\n")
cat("Observed max span:", obs_max, "windows (", obs_max * 0.1, "Mb)\n")
cat("Null distribution summary:\n")
cat("  Min:", min(null_max_spans), "\n")
cat("  Median:", median(null_max_spans), "\n")
cat("  95th percentile:", quantile(null_max_spans, 0.95), "\n")
cat("  99th percentile:", quantile(null_max_spans, 0.99), "\n")
cat("  Max:", max(null_max_spans), "\n")
cat("P-value:", ifelse(p_value == 0, paste0("< ", 1/n_perm), p_value), "\n")

# =============================================================================
# Final results table
# =============================================================================

final_table <- results %>%
  mutate(
    second_sample = str_extract(pair, "(?<=:).*"),
    farm = str_extract(second_sample, "(?<=Treatment)\\d+"),
    timepoint = str_extract(second_sample, "(Pre|Post)$"),
    Comparison = paste0("Farm ", farm, " ", timepoint),
    `Block Coordinates (Mb)` = paste0(start_mb, " – ", end_mb),
    `Width (Mb)` = width_mb
  ) %>%
  select(Comparison, `Block Coordinates (Mb)`, `Width (Mb)`, 
         `Span (windows)` = span_windows, `Windows > threshold` = n_pass) %>%
  arrange(
    as.numeric(str_extract(Comparison, "\\d+")),
    desc(str_detect(Comparison, "Pre"))
  )

cat("\n=== Final Results Table ===\n")
cat("Threshold: 99th percentile FST =", round(cutoff, 4), "\n")
cat("Beta-tubulin position: 7.025 Mb\n\n")
print(final_table)

# Save table
write_csv(final_table, "selection_blocks_results.csv")
cat("\nTable saved to: selection_blocks_results.csv\n")

# =============================================================================
# Save permutation results
# =============================================================================

perm_summary <- tibble(
  parameter = c("threshold", "observed_max_span", "null_median", "null_95th", 
                "null_99th", "null_max", "p_value", "n_permutations"),
  value = c(round(cutoff, 4), obs_max, median(null_max_spans), 
            quantile(null_max_spans, 0.95), quantile(null_max_spans, 0.99),
            max(null_max_spans), p_value, n_perm)
)

write_csv(perm_summary, "permutation_test_results.csv")
cat("Permutation results saved to: permutation_test_results.csv\n")

# Clean up parallel workers
plan(sequential)



# Compare FST distributions: autosomes vs X
fst_100kb %>%
  filter(pair %in% signal_pairs) %>%
  mutate(
    chrom_type = case_when(
      str_detect(chrom, "chrX") ~ "X",
      str_detect(chrom, "chr[1-5]") ~ "Autosome",
      TRUE ~ "Other"
    )
  ) %>%
  filter(chrom_type %in% c("Autosome", "X")) %>%
  group_by(chrom_type) %>%
  summarise(
    mean_fst = mean(fst, na.rm = TRUE),
    median_fst = median(fst, na.rm = TRUE),
    q99 = quantile(fst, 0.99, na.rm = TRUE),
    n = n()
  )





# Threshold with X included
bg_data_with_X <- fst_100kb %>%
  filter(pair %in% signal_pairs) %>%
  filter(str_detect(chrom, "chr[1-5X]"))

cutoff_with_X <- quantile(bg_data_with_X$fst, 0.99, na.rm = TRUE)

cat("Threshold without X:", round(cutoff, 4), "\n")
cat("Threshold with X:", round(cutoff_with_X, 4), "\n")

# Check for blocks on X
x_chrom <- "hcontortus_chrX_Celeg_TT_arrow_pilon"

x_blocks <- map_dfr(signal_pairs, function(p) {
  blocks <- fst_100kb %>%
    filter(pair == p, chrom == x_chrom) %>%
    find_all_blocks(cutoff, max_gap = 3) %>%
    filter(span_windows >= 10)
  
  if (nrow(blocks) > 0) {
    blocks %>% mutate(pair = p)
  } else {
    tibble()
  }
})

cat("\nBlocks on X chromosome:\n")
if (nrow(x_blocks) > 0) {
  print(x_blocks)
} else {
  cat("None\n")
}







# Compare blocks with both thresholds
cat("=== Blocks with autosome-only threshold (0.1139) ===\n")
print(results)

cat("\n=== Blocks with X-included threshold (0.1174) ===\n")
results_with_X <- map_dfr(signal_pairs, function(p) {
  blocks <- fst_100kb %>%
    filter(pair == p, chrom == target_chrom) %>%
    find_all_blocks(cutoff_with_X, max_gap = 3) %>%
    filter(span_windows >= 10)
  
  if (nrow(blocks) > 0) {
    blocks %>% 
      mutate(pair = p, width_mb = round(end_mb - start_mb, 2)) %>%
      select(pair, start_mb, end_mb, width_mb, span_windows, n_pass)
  } else {
    tibble(pair = p, start_mb = NA, end_mb = NA, width_mb = NA, 
           span_windows = NA, n_pass = NA)
  }
})

print(results_with_X)

# Side-by-side comparison
cat("\n=== Comparison ===\n")
comparison <- results %>%
  select(pair, span_autosome = span_windows, width_autosome = width_mb) %>%
  left_join(
    results_with_X %>% select(pair, span_with_X = span_windows, width_with_X = width_mb),
    by = "pair"
  )
print(comparison)





# =============================================================================
# Plot with X-included threshold (0.1174)
# =============================================================================

cutoff <- cutoff_with_X  # 0.1174

# Get all pairs starting with Treatment02_Pre
all_pairs <- fst_100kb %>%
  filter(chrom == target_chrom) %>%
  pull(pair) %>%
  unique() %>%
  str_subset("^Treatment02_Pre:")

# Get blocks with new threshold
all_blocks <- map_dfr(all_pairs, function(p) {
  blocks <- fst_100kb %>%
    filter(pair == p, chrom == target_chrom) %>%
    find_all_blocks(cutoff, max_gap = 3) %>%
    filter(span_windows >= 10)
  
  if (nrow(blocks) > 0) {
    blocks %>% mutate(pair = p)
  } else {
    tibble()
  }
})

# Prepare FST data
fst_chr1_all <- fst_100kb %>%
  filter(chrom == target_chrom, pair %in% signal_pairs) %>%
  mutate(pos_mb = midpoint / 1e6)

# Add in_block flag
if (nrow(all_blocks) > 0) {
  fst_chr1_all <- fst_chr1_all %>%
    left_join(
      all_blocks %>% select(pair, start_mb, end_mb),
      by = "pair",
      relationship = "many-to-many"
    ) %>%
    mutate(in_block = pos_mb >= start_mb & pos_mb <= end_mb) %>%
    group_by(pair, pos_mb, fst) %>%
    summarise(in_block = any(in_block, na.rm = TRUE), .groups = "drop") %>%
    mutate(in_block = replace_na(in_block, FALSE))
} else {
  fst_chr1_all <- fst_chr1_all %>%
    mutate(in_block = FALSE)
}

# Create labels
fst_chr1_all <- fst_chr1_all %>%
  mutate(
    second_sample = str_extract(pair, "(?<=:).*"),
    farm_num = as.numeric(str_extract(second_sample, "(?<=Treatment)\\d+")),
    timepoint = str_extract(second_sample, "(Pre|Post)$"),
    pair_label = paste0("Farm ", sprintf("%02d", farm_num), " ", timepoint)
  )

# Order by farm, Pre before Post
farm_order <- fst_chr1_all %>%
  distinct(farm_num, timepoint, pair_label) %>%
  arrange(farm_num, desc(timepoint == "Pre")) %>%
  pull(pair_label)

fst_chr1_all <- fst_chr1_all %>%
  mutate(pair_label = factor(pair_label, levels = farm_order))

# Plot
p_chr1 <- ggplot(fst_chr1_all, aes(x = pos_mb, y = fst, color = in_block)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "gray50", "TRUE" = "firebrick"),
    labels = c("FALSE" = "Outside block", "TRUE" = "In significant block"),
    name = NULL
  ) +
  geom_vline(xintercept = beta_tubulin_pos / 1e6, color = "darkgreen", 
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dotted", linewidth = 0.5) +
  facet_wrap(~ pair_label, ncol = 1) +
  labs(
    x = "Position on chromosome I (Mb)",
    y = expression(F[ST]),
    title = "Selection blocks on chromosome I",
    subtitle = paste0("Threshold: 99th percentile FST = ", round(cutoff, 4), 
                      " (autosomes + X)")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 8)
  )

print(p_chr1)

# Save
ggsave("chr1_selection_blocks_with_X.pdf", p_chr1, width = 10, height = 12)
cat("Plot saved to: chr1_selection_blocks_with_X.pdf\n")








# =============================================================================
# Genome-wide plot (Farm 11 Pre as example)
# =============================================================================

focus_pair <- "Treatment02_Pre:Treatment11_Pre"

# Get blocks across all chromosomes
all_chroms <- fst_100kb %>%
  filter(!str_detect(chrom, "mito")) %>%
  pull(chrom) %>%
  unique()

genome_blocks <- map_dfr(all_chroms, function(chr) {
  blocks <- fst_100kb %>%
    filter(pair == focus_pair, chrom == chr) %>%
    find_all_blocks(cutoff, max_gap = 3) %>%
    filter(span_windows >= 10)
  
  if (nrow(blocks) > 0) {
    blocks %>% mutate(chrom = chr)
  } else {
    tibble()
  }
})

cat("Blocks found genome-wide (Farm 11 Pre):\n")
if (nrow(genome_blocks) > 0) {
  genome_blocks %>%
    mutate(chrom_short = str_extract(chrom, "chr[0-9X]+")) %>%
    select(chrom_short, start_mb, end_mb, span_windows, n_pass) %>%
    print()
} else {
  cat("None\n")
}

# Prepare genome data
fst_genome <- fst_100kb %>%
  filter(pair == focus_pair, !str_detect(chrom, "mito")) %>%
  mutate(pos_mb = midpoint / 1e6)

# Add in_block flag
if (nrow(genome_blocks) > 0) {
  fst_genome <- fst_genome %>%
    left_join(
      genome_blocks %>% select(chrom, start_mb, end_mb),
      by = "chrom",
      relationship = "many-to-many"
    ) %>%
    mutate(in_block = pos_mb >= start_mb & pos_mb <= end_mb) %>%
    group_by(chrom, pos_mb, fst) %>%
    summarise(in_block = any(in_block, na.rm = TRUE), .groups = "drop") %>%
    mutate(in_block = replace_na(in_block, FALSE))
} else {
  fst_genome <- fst_genome %>%
    mutate(in_block = FALSE)
}

# Create chromosome labels
fst_genome <- fst_genome %>%
  mutate(
    chrom_num = str_extract(chrom, "chr(\\d+|X)"),
    chrom_num = str_replace(chrom_num, "chr", "Chr ")
  )

chrom_order <- c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr X")
fst_genome <- fst_genome %>%
  mutate(chrom_num = factor(chrom_num, levels = chrom_order))

# Beta-tubulin marker (chr 1 only)
beta_tub_df <- tibble(
  chrom_num = factor("Chr 1", levels = chrom_order),
  x = beta_tubulin_pos / 1e6
)

# Plot
p_genome <- ggplot(fst_genome, aes(x = pos_mb, y = fst, color = in_block)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(
    values = c("FALSE" = "gray50", "TRUE" = "firebrick"),
    labels = c("FALSE" = "Outside block", "TRUE" = "In significant block"),
    name = NULL
  ) +
  geom_vline(data = beta_tub_df, aes(xintercept = x),
             color = "darkgreen", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dotted", linewidth = 0.5) +
  facet_wrap(~ chrom_num, ncol = 1, scales = "free_x") +
  labs(
    x = "Position (Mb)",
    y = expression(F[ST]),
    title = "Genome-wide FST: Farm 11 Pre vs Baseline",
    subtitle = paste0("Threshold: 99th percentile = ", round(cutoff, 4))
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 9)
  )

print(p_genome)
ggsave("genome_wide_fst_farm11pre.pdf", p_genome, width = 10, height = 10)

# =============================================================================
# Final results table with X-included threshold
# =============================================================================

final_table <- results_with_X %>%
  mutate(
    second_sample = str_extract(pair, "(?<=:).*"),
    farm = str_extract(second_sample, "(?<=Treatment)\\d+"),
    timepoint = str_extract(second_sample, "(Pre|Post)$"),
    Comparison = paste0("Farm ", farm, " ", timepoint)
  ) %>%
  select(
    Comparison,
    `Start (Mb)` = start_mb,
    `End (Mb)` = end_mb,
    `Width (Mb)` = width_mb,
    `Windows` = span_windows,
    `Windows > threshold` = n_pass
  ) %>%
  arrange(
    as.numeric(str_extract(Comparison, "\\d+")),
    desc(str_detect(Comparison, "Pre"))
  )

cat("\n=== Final Results Table ===\n")
cat("Threshold: 99th percentile FST =", round(cutoff, 4), "(autosomes + X)\n")
cat("Beta-tubulin position: 7.025 Mb\n\n")
print(final_table)

write_csv(final_table, "selection_blocks_final.csv")
cat("\nSaved to: selection_blocks_final.csv\n")