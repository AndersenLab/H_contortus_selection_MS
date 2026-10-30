#setwd(dir = "~/Etta/H_contortus_selection/scripts/")

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

# --------------------------
# --------------------------
# 1. Data loading functions 
# --------------------------
# --------------------------

# --------------------------
# Define load_diversity()
# Input: Raw diversity output from grenedalf (diversity.csv)
# Output: A list of three dataframes
# 1. windows: The number of sites that were included or excluded by window
# 2. qc: The number of sites that were included or excluded by window and by sample
# 3. values: The numeric output (pi, theta, tajimas_d) by window and by sample

load_diversity <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE) # Read raw data and suppress warning messages
  
  # Select columns for windows dataframe
  windows <- raw %>%
    dplyr::select(chrom, start, end, starts_with("total."))
  
  # Each sample is a column prefix (Example: Treatment09_Pre.theta_pi)
  # Convert to long format where each sample gets its own row
  long <- raw %>%
    dplyr::select(chrom, start, end, matches("^[^.]+\\.(missing|numeric|passed|theta_pi|theta_watterson|tajimas_d)$")) %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = c("sample", "metric"),
      names_pattern = "^(.+)\\.(missing|numeric|passed|theta_pi|theta_watterson|tajimas_d)$",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value)) # Convert NaN to NA
  
  # Pivot metrics back out into their own columns, one row per window per sample
  wide <- long %>%
    pivot_wider(names_from = metric, values_from = value)
  
  # Make qc dataframe
  qc <- wide %>%
    dplyr::select(chrom, start, end, sample, missing, numeric, passed)
  
  # Make diversity values dataframe
  values <- wide %>%
    dplyr::select(chrom, start, end, sample, theta_pi, theta_watterson, tajimas_d)
  
  # Store as a list
  list(
    windows = windows,
    qc = qc,
    values = values
  )
}
# --------------------------

# --------------------------
# --------------------------
# 2. Data processing functions 
# --------------------------
# --------------------------

# --------------------------
# Define relabel_sample()
# Converts grenedalf naming structure to Pre_X or Post_X for plotting

relabel_sample <- function(x) {
  timepoint <- ifelse(grepl("_Pre$", x), "Pre", "Post")
  number <- as.integer(str_extract(x, "(?<=Treatment)[0-9]+"))
  paste0(timepoint, "_", number)
}
# --------------------------

# --------------------------
# Define normalize_and_calculate_pi_ratio()
# Input: Diversity values from load_diversity(), timepoint to use as baseline
# Output: Pi ratio (compared to baseline timepoint) for each subsequent timepoint, both raw and log10. Also the genome-wide rank for that timepoint.

normalize_and_calculate_pi_ratio <- function(data, baseline = "F0") {
  
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

# --------------------------
# Define identify_pi_outliers()
# Input: data (output of calculate_pi_ratio()), percentile cutoff for "low diversity"
# Output: same data, with an outlier flag column added, using one threshold shared across all samples

identify_pi_outliers <- function(data, percentile = 0.01) {
  threshold <- quantile(data$pi_ratio_log10, probs = percentile, na.rm = TRUE)
  
  data %>%
    mutate(
      threshold = threshold,
      is_outlier = pi_ratio_log10 <= threshold
    )
}
# --------------------------

# --------------------------
# Define find_contiguous_regions()
# Input: data (output of identify_pi_outliers()), minimum consecutive outlier windows to call a region
# Output: one row per candidate region, per sample/chromosome, with its span and window count

find_contiguous_regions <- function(data, min_windows = 3) {
  
  data %>%
    arrange(sample, chrom, start) %>%
    group_by(sample, chrom) %>%
    mutate(
      run_id = cumsum(is_outlier != lag(is_outlier, default = first(is_outlier)))
    ) %>%
    dplyr::filter(is_outlier) %>%
    group_by(sample, chrom, run_id) %>%
    summarize(
      region_start = min(start),
      region_end = max(end),
      n_windows = n(),
      mean_pi_ratio_log10 = mean(pi_ratio_log10, na.rm = TRUE),
      min_pi_ratio_log10 = min(pi_ratio_log10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_windows >= min_windows) %>%
    dplyr::select(-run_id) %>%
    arrange(sample, chrom, region_start)
}
# --------------------------

# --------------------------
# Define flag_region_windows()
# Input: data (output of identify_pi_outliers()), candidate regions (output of find_contiguous_regions())
# Output: same data, with a new column marking whether each window falls inside a validated candidate region

flag_region_windows <- function(data, regions) {
  data %>%
    left_join(
      regions %>% dplyr::select(sample, chrom, region_start, region_end),
      by = c("sample", "chrom"),
      relationship = "many-to-many"
    ) %>%
    group_by(sample, chrom, start, end) %>%
    mutate(in_region = any(start >= region_start & end <= region_end, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(sample, chrom, start, end, .keep_all = TRUE) %>%
    mutate(in_region = tidyr::replace_na(in_region, FALSE))
}
# --------------------------

# --------------------------
# --------------------------
# 3. Plotting functions 
# --------------------------
# --------------------------

# --------------------------
# Define plot_diversity()
# Input: grenedalf dataframe with all three metrics, which metric you want to plot, and which timepoints you want to plot
# Output: Genome-wide plot, faceted by sample, and colored by chromosome

plot_diversity <- function(data,
                           samples,
                           metric = "theta_pi",
                           chrom_order = c(
                             "hcontortus_chr1_Celeg_TT_arrow_pilon",
                             "hcontortus_chr2_Celeg_TT_arrow_pilon",
                             "hcontortus_chr3_Celeg_TT_arrow_pilon",
                             "hcontortus_chr4_Celeg_TT_arrow_pilon",
                             "hcontortus_chr5_Celeg_TT_arrow_pilon",
                             "hcontortus_chrX_Celeg_TT_arrow_pilon"
                           ),
                           short_labels = c("I", "II", "III", "IV", "V", "X"),
                           show_label = FALSE,
                           smooth = FALSE,
                           loess_span = 0.1,
                           save_path = NULL,
                           width = 10,
                           height = 3,
                           dpi = 300) {
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(.data[[metric]])) %>%
    mutate(chrom = factor(chrom, levels = chrom_order))
  
  # Build cumulative genome-wide offset for each chromosome, so positions
  # concatenate end-to-end rather than each starting back at zero
  chrom_lengths <- plot_data %>%
    group_by(chrom) %>%
    summarize(max_end = max(end), .groups = "drop") %>%
    arrange(match(chrom, chrom_order)) %>%
    mutate(offset = lag(cumsum(as.numeric(max_end)), default = 0))
  
  plot_data <- plot_data %>%
    left_join(chrom_lengths %>% dplyr::select(chrom, offset), by = "chrom") %>%
    mutate(
      position_mb = (((start + end) / 2) + offset) / 1e6,
      chrom_short = factor(short_labels[match(chrom, chrom_order)], levels = short_labels)
    )
  
  # Chromosome midpoints, for x-axis tick labels
  axis_labels <- plot_data %>%
    group_by(chrom_short) %>%
    summarize(mid = mean(range(position_mb)), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = position_mb, y = .data[[metric]], color = chrom_short)) +
    geom_point(alpha = 0.5, size = 0.6) +
    scale_color_viridis_d() +
    scale_x_continuous(breaks = axis_labels$mid, labels = axis_labels$chrom_short) +
    labs(x = "Chromosome", y = expression(Nucleotide~diversity~(pi)), color = "Chromosome") +
    guides(color = "none") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  # Add facet with or without labels
  if (show_label) {
    p <- p + facet_wrap(~sample, ncol = 1, labeller = as_labeller(relabel_sample))
  } else {
    p <- p + 
      facet_wrap(~sample, ncol = 1) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank()
      )
  }
  
  if (smooth) {
    p <- p + geom_smooth(aes(group = chrom_short), method = "loess", span = loess_span,
                         se = FALSE, color = "black", linewidth = 0.5)
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height * length(samples), dpi = dpi)
  }
  
  p
}

# --------------------------

# --------------------------
# Define plot_pi_ratio()
# Input: data (output of calculate_pi_ratio()), one or more samples to plot
# Output: Genome-wide plot of log10(pi ratio) to baseline, faceted by chromosome and sample
plot_pi_ratio <- function(data,
                          samples,
                          compact = FALSE,
                          roman_numerals = TRUE,
                          chrom_order = c(
                            "hcontortus_chr1_Celeg_TT_arrow_pilon",
                            "hcontortus_chr2_Celeg_TT_arrow_pilon",
                            "hcontortus_chr3_Celeg_TT_arrow_pilon",
                            "hcontortus_chr4_Celeg_TT_arrow_pilon",
                            "hcontortus_chr5_Celeg_TT_arrow_pilon",
                            "hcontortus_chrX_Celeg_TT_arrow_pilon"
                          ),
                          save_path = NULL,
                          width = 7.5,
                          height_per_sample = 2,
                          dpi = 300) {
  
  # Set chromosome labels based on roman_numerals
  if (roman_numerals) {
    short_labels <- c("I", "II", "III", "IV", "V", "X")
  } else {
    short_labels <- c("1", "2", "3", "4", "5", "X")
  }
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(pi_ratio_log10)) %>%
    mutate(
      chrom = factor(chrom, levels = chrom_order),
      chrom_short = factor(short_labels[match(chrom, chrom_order)], levels = short_labels),
      position_mb = ((start + end) / 2) / 1e6
    )
  
  # Relabel samples to human-readable names
  if (compact) {
    plot_data <- plot_data %>%
      mutate(sample_label = gsub("Treatment0?(\\d+)_(Pre|Post)", "\\2-Txt \\1", sample),
             sample_label = factor(sample_label, levels = unique(sample_label)))
  } else {
    plot_data <- plot_data %>%
      mutate(sample_label = gsub("Treatment0?(\\d+)_(Pre|Post)", "\\2-treatment \\1", sample),
             sample_label = factor(sample_label, levels = unique(sample_label)))
  }
  
  # Standard theme
  my_theme <- theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 0)),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10, margin = margin(r = 0)),
      strip.text = element_text(size = 12),
      strip.text.y.left = element_text(angle = 0, size = 12, margin = margin(l = 0, r = 0)),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.1, "lines")
    )
  
  # Set x-axis breaks
  if (compact) {
    x_scale <- scale_x_continuous(breaks = function(x) setdiff(pretty(x), 0), expand = expansion(mult = c(0.02, 0.02)))
  } else {
    x_scale <- scale_x_continuous(breaks = c(15, 30, 45), expand = expansion(mult = c(0.02, 0.02)))
  }
  
  p <- ggplot(plot_data, aes(x = position_mb, y = pi_ratio_log10)) +
    geom_point(alpha = 0.3, size = 0.6, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    x_scale +
    labs(x = "Position (Mb)", y = "**log<sub>10</sub> (π ratio)**") +
    my_theme +
    theme(axis.title.y = ggtext::element_markdown(size = 12))
  
  # Facet by chromosome (columns) and sample (rows)
  if (length(samples) == 1) {
    p <- p + 
      facet_wrap(~chrom_short, nrow = 1, scales = "free_x")
  } else {
    p <- p + 
      facet_grid(sample_label ~ chrom_short, scales = "free_x", switch = "y") +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 12)
      )
  }
  
  if (!is.null(save_path)) {
    total_height <- height_per_sample * length(samples)
    ggsave(save_path, plot = p, width = width, height = total_height, dpi = dpi)
  }
  
  p
}

# --------------------------
# Define plot_pi_ratio_by_value()
# Input: data (output of calculate_pi_ratio()), one or more samples to plot
# Output: Genome-wide plot colored by pi ratio value instead of chromosome

plot_pi_ratio_by_value <- function(data,
                                   samples,
                                   chrom_order,
                                   save_path = NULL,
                                   width = 10,
                                   height = 3,
                                   dpi = 300) {
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(pi_ratio_log10)) %>%
    mutate(position_mb = (start + end) / 2 / 1e6)
  
  p <- ggplot(plot_data, aes(x = position_mb, y = pi_ratio_log10, color = pi_ratio_log10)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_viridis_c() +
    facet_wrap(~sample, ncol = 1, strip.position = "left", labeller = as_labeller(relabel_sample)) +
    labs(x = "Position (Mb)", y = expression(log[10]~"("~pi~"ratio)"), color = expression(log[10]~"ratio")) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0)
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height * length(samples), dpi = dpi)
  }
  
  p
}
# --------------------------

# --------------------------
# Define plot_pi_outliers()
# Input: data (output of flag_region_windows()), one or more samples to plot
# Output: Genome-wide Manhattan-style plot of log10(pi ratio) to baseline, faceted by sample, colored by validated region membership

plot_pi_outliers <- function(data,
                             samples,
                             chrom_order = c(
                               "hcontortus_chr1_Celeg_TT_arrow_pilon",
                               "hcontortus_chr2_Celeg_TT_arrow_pilon",
                               "hcontortus_chr3_Celeg_TT_arrow_pilon",
                               "hcontortus_chr4_Celeg_TT_arrow_pilon",
                               "hcontortus_chr5_Celeg_TT_arrow_pilon",
                               "hcontortus_chrX_Celeg_TT_arrow_pilon"
                             ),
                             short_labels = c("I", "II", "III", "IV", "V", "X"),
                             save_path = NULL,
                             width = 10,
                             height = 3,
                             dpi = 300) {
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(pi_ratio_log10)) %>%
    mutate(chrom = factor(chrom, levels = chrom_order))
  
  chrom_lengths <- plot_data %>%
    group_by(chrom) %>%
    summarize(max_end = max(end), .groups = "drop") %>%
    arrange(match(chrom, chrom_order)) %>%
    mutate(offset = lag(cumsum(as.numeric(max_end)), default = 0))
  
  plot_data <- plot_data %>%
    left_join(chrom_lengths %>% dplyr::select(chrom, offset), by = "chrom") %>%
    mutate(
      position_mb = (((start + end) / 2) + offset) / 1e6,
      chrom_short = factor(short_labels[match(chrom, chrom_order)], levels = short_labels)
    )
  
  axis_labels <- plot_data %>%
    group_by(chrom_short) %>%
    summarize(mid = mean(range(position_mb)), .groups = "drop")
  
  # Plot background windows first, validated-region windows last (drawn on top)
  p <- ggplot(plot_data %>% arrange(in_region),
              aes(x = position_mb, y = pi_ratio_log10, color = in_region)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "firebrick"),
                       labels = c("FALSE" = "Background", "TRUE" = "Candidate region")) +
    scale_x_continuous(breaks = axis_labels$mid, labels = axis_labels$chrom_short) +
    facet_wrap(~sample, ncol = 1, strip.position = "left", labeller = as_labeller(relabel_sample)) +
    labs(x = "Chromosome", y = expression(log[10]~"("~pi~"ratio)"), color = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0)
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height * length(samples), dpi = dpi)
  }
  
  p
}
# --------------------------

# --------------------------
# --------------------------
# 4. Analysis pipeline 
# --------------------------
# --------------------------

# --------------------------
# Load data
diversity_data <- load_diversity("~/Etta/H_contortus_selection/diversity/10X/diversity.csv")

# View dataframes
diversity_data$windows
diversity_data$qc
diversity_data$values

# Save values to an object
diversity_values <- diversity_data$values

# Check available samples
unique(diversity_values$sample)
# --------------------------

# --------------------------
# Calculate pi ratios
pi_ratios <- normalize_and_calculate_pi_ratio(diversity_values, baseline = "Treatment02_Pre")

head(pi_ratios)

# Get all sample names for plotting
all_samples <- unique(pi_ratios$sample)
# --------------------------

# --------------------------
# Identify outliers and candidate regions
pi_outliers <- identify_pi_outliers(pi_ratios, percentile = 0.001)

candidate_regions <- find_contiguous_regions(pi_outliers, min_windows = 3)
candidate_regions

pi_flagged <- flag_region_windows(pi_outliers, candidate_regions)
# --------------------------

# --------------------------
# --------------------------
# 5. Generate plots 
# --------------------------
# --------------------------

# Define output directory
out_dir <- "/vast/eande106/projects/Etta/H_contortus_selection/figures/"

# --------------------------
# Plot raw diversity
plot_diversity(
  diversity_values,
  samples = "Treatment11_Pre",
  metric = "theta_pi"
)
ggsave(paste0(out_dir, "diversity_Treatment11_Pre.jpg"), plot = last_plot(), width = 3.75, height = 2.5, dpi = 300)
# --------------------------

# --------------------------
# Plot pi ratio for final sample
plot_pi_ratio(
  pi_ratios,
  samples = "Treatment11_Pre"
)
ggsave(paste0(out_dir, "pi_ratio_Treatment11_Pre.jpg"), plot = last_plot(), width = 3.75, height = 2.5, dpi = 300)
# --------------------------

# --------------------------
# Plot pi ratio for all samples

# Full labels, roman numerals (default)
plot_pi_ratio(
  pi_ratios,
  samples = c("Treatment03_Pre", "Treatment04_Pre", "Treatment08_Pre",
              "Treatment08_Post", "Treatment09_Pre", "Treatment09_Post",
              "Treatment10_Pre", "Treatment10_Post", "Treatment11_Pre")
)
ggsave(paste0(out_dir, "pi_ratio_roman.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)


# Compact labels, roman numerals
plot_pi_ratio(
  pi_ratios,
  samples = c("Treatment03_Pre", "Treatment04_Pre", "Treatment08_Pre",
              "Treatment08_Post", "Treatment09_Pre", "Treatment09_Post",
              "Treatment10_Pre", "Treatment10_Post", "Treatment11_Pre"),
  compact = TRUE
)
ggsave(paste0(out_dir, "pi_ratio_compact_roman.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)

# Full labels, arabic numerals
plot_pi_ratio(
  pi_ratios,
  samples = c("Treatment03_Pre", "Treatment04_Pre", "Treatment08_Pre",
              "Treatment08_Post", "Treatment09_Pre", "Treatment09_Post",
              "Treatment10_Pre", "Treatment10_Post", "Treatment11_Pre"),
  roman_numerals = FALSE
)
ggsave(paste0(out_dir, "pi_ratio_arabic.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)


# Compact labels, arabic numerals
plot_pi_ratio(
  pi_ratios,
  samples = c("Treatment03_Pre", "Treatment04_Pre", "Treatment08_Pre",
              "Treatment08_Post", "Treatment09_Pre", "Treatment09_Post",
              "Treatment10_Pre", "Treatment10_Post", "Treatment11_Pre"),
  compact = TRUE,
  roman_numerals = FALSE
)
ggsave(paste0(out_dir, "pi_ratio_compact_arabic.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)
# --------------------------

# --------------------------
# Plot pi ratio colored by value (single chromosome)
plot_pi_ratio_by_value(
  pi_ratios,
  samples = all_samples,
  chrom_order = "hcontortus_chr1_Celeg_TT_arrow_pilon"
)
ggsave(paste0(out_dir, "pi_ratio_by_value_chr1.jpg"), plot = last_plot(), width = 3.75, height = 8, dpi = 300)
# --------------------------

# --------------------------
# Plot outliers
plot_pi_outliers(
  pi_flagged,
  samples = all_samples
)
ggsave(paste0(out_dir, "pi_outliers_all_samples.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)

