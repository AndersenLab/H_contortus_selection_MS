# --------------------------
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtext)
library(patchwork)
library(viridisLite)
# --------------------------

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

# --------------------------
# Define load_fst()
# Input: Raw FST output from grenedalf (fst.csv)
# Output: A list of three dataframes (windows, qc, values)

load_fst <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  
  windows <- raw %>%
    dplyr::select(chrom, start, end, starts_with("total."))
  
  long <- raw %>%
    dplyr::select(chrom, start, end, matches("^[^.]+:[^.]+\\.(missing|numeric|passed|fst)$")) %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = c("sample_pair", "metric"),
      names_pattern = "^(.+)\\.(missing|numeric|passed|fst)$",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value))
  
  wide <- long %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    separate(sample_pair, into = c("sample_1", "sample_2"), sep = ":", remove = FALSE)
  
  qc <- wide %>%
    dplyr::select(chrom, start, end, sample_pair, sample_1, sample_2, missing, numeric, passed)
  
  values <- wide %>%
    dplyr::select(chrom, start, end, sample_pair, sample_1, sample_2, fst)
  
  list(windows = windows, qc = qc, values = values)
}
# --------------------------

# --------------------------
# Define normalize_and_calculate_pi_ratio()
# Input: Diversity values from load_diversity(), timepoint to use as baseline
# Output: Pi ratio (compared to baseline timepoint) for each subsequent timepoint

normalize_and_calculate_pi_ratio <- function(data, baseline = "F0") {
  
  normalized_data <- data %>%
    group_by(sample) %>%
    mutate(pi_norm = theta_pi / median(theta_pi, na.rm = TRUE)) %>%
    ungroup()
  
  baseline_data <- normalized_data %>%
    dplyr::filter(sample == baseline) %>%
    dplyr::select(chrom, start, end, pi_norm_baseline = pi_norm)
  
  other_data <- normalized_data %>%
    dplyr::filter(sample != baseline) %>%
    dplyr::select(chrom, start, end, sample, pi_norm_sample = pi_norm)
  
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
# --- Comparison 1: both at 10kb ---
diversity_data_10kb <- load_diversity("~/Etta/H_contortus_selection/diversity/genome_10kb/diversity.csv")
diversity_values_10kb <- diversity_data_10kb$values
pi_ratios_10kb <- normalize_and_calculate_pi_ratio(diversity_values_10kb, baseline = "Treatment02_Pre")

fst_data_10kb <- load_fst("~/Etta/H_contortus_selection/fst/10X/fst.csv")  # original 10kb FST run
fst_values_10kb <- fst_data_10kb$values

# --- Comparison 2: both at 100kb ---
diversity_data_100kb <- load_diversity("~/Etta/H_contortus_selection/diversity/10X/diversity.csv")  # original 100kb diversity run
diversity_values_100kb <- diversity_data_100kb$values
pi_ratios_100kb <- normalize_and_calculate_pi_ratio(diversity_values_100kb, baseline = "Treatment02_Pre")

fst_data_100kb <- load_fst("~/Etta/H_contortus_selection/fst/genome_100kb/fst.csv")  # new 100kb FST run
fst_values_100kb <- fst_data_100kb$values
# --------------------------

# --------------------------
# Reference region (isotype-1 beta-tubulin)
tbb_iso1_region <- c(7.027492, 7.031447)
target_chrom <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
# --------------------------

# --------------------------
# Define plot_pi_tajima_fst_combined()
# Input: pi ratio values, fst values, pi-ratio sample name, sample_pair name (fst), chromosome
# Optional: diversity values and tajima sample for Tajima's D panel
# Output: stacked pi ratio + (Tajima's D) + FST plot, shared x-axis, gene region highlighted

plot_pi_tajima_fst_combined <- function(pi_ratio_values,
                                        fst_values,
                                        pi_sample,
                                        fst_pair,
                                        chrom,
                                        diversity_values = NULL,
                                        tajima_sample = NULL,
                                        show_tajima = FALSE,
                                        gene_region = NULL,
                                        xlim_range = NULL,
                                        loess_span = 0.1,
                                        roman_numerals = FALSE,
                                        save_path = NULL,
                                        width = 7.5,
                                        height = 5,
                                        dpi = 300) {
  
  # Viridis blue color for points
  point_color <- viridisLite::viridis(10)[4]
  
  # Shared theme
  my_theme <- theme_bw(base_size = 12) +
    theme(
      axis.title.x = ggtext::element_markdown(face = "bold", size = 12),
      axis.title.y = ggtext::element_markdown(face = "bold", size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
  
  # Gene midpoint for asterisk
  gene_mid <- if (!is.null(gene_region)) mean(gene_region) else NULL
  
  # --- Pi ratio panel ---
  pi_data <- pi_ratio_values %>%
    dplyr::filter(sample == pi_sample, chrom == !!chrom, !is.na(pi_ratio_log10)) %>%
    mutate(position_mb = ((start + end) / 2) / 1e6)
  
  p_pi <- ggplot(pi_data, aes(x = position_mb, y = pi_ratio_log10)) +
    geom_point(alpha = 0.4, size = 0.6, color = point_color) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = NULL, y = "**log<sub>10</sub>(π ratio)**") +
    my_theme
  
  if (!is.null(gene_mid)) {
    p_pi <- p_pi + annotate("text", x = gene_mid, y = Inf, label = "*", 
                            vjust = 1.5, size = 8, color = "black")
  }
  
  # --- Tajima's D panel (optional) ---
  p_taj <- NULL
  if (show_tajima) {
    if (is.null(diversity_values) || is.null(tajima_sample)) {
      stop("diversity_values and tajima_sample required when show_tajima = TRUE")
    }
    
    taj_data <- diversity_values %>%
      dplyr::filter(sample == tajima_sample, chrom == !!chrom, !is.na(tajimas_d)) %>%
      mutate(position_mb = ((start + end) / 2) / 1e6)
    
    p_taj <- ggplot(taj_data, aes(x = position_mb, y = tajimas_d)) +
      geom_point(alpha = 0.4, size = 0.6, color = point_color) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      labs(x = NULL, y = "**Tajima's D**") +
      my_theme
    
    if (!is.null(gene_mid)) {
      p_taj <- p_taj + annotate("text", x = gene_mid, y = Inf, label = "*", 
                                vjust = 1.5, size = 8, color = "black")
    }
  }
  
  # --- FST panel ---
  fst_data_chr <- fst_values %>%
    dplyr::filter(sample_pair == fst_pair, chrom == !!chrom, !is.na(fst)) %>%
    mutate(position_mb = ((start + end) / 2) / 1e6)
  
  p_fst <- ggplot(fst_data_chr, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.4, size = 0.6, color = point_color) +
    geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = "**Position (Mb)**", y = "**F<sub>ST</sub>**") +
    my_theme
  
  if (!is.null(gene_mid)) {
    p_fst <- p_fst + annotate("text", x = gene_mid, y = Inf, label = "*", 
                              vjust = 1.5, size = 8, color = "black")
  }
  
  # --- Combine panels ---
  if (show_tajima) {
    combined <- p_pi / p_taj / p_fst + plot_layout(axes = "collect")
  } else {
    combined <- p_pi / p_fst + plot_layout(axes = "collect")
  }
  
  if (!is.null(xlim_range)) {
    combined <- combined & xlim(xlim_range[1], xlim_range[2])
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = combined, width = width, height = height, dpi = dpi)
  }
  
  combined
}
# --------------------------

# --------------------------
out_dir <- "~/Etta/H_contortus_selection/figures/"
# --------------------------

p_10kb <- plot_pi_tajima_fst_combined(
  pi_ratios_10kb,
  fst_values_10kb,
  pi_sample = "Treatment11_Pre",
  fst_pair = "Treatment02_Pre:Treatment11_Pre",
  chrom = target_chrom,
  gene_region = tbb_iso1_region  
) + plot_annotation(tag_levels = 'A')


p_10kb
ggsave(paste0(out_dir, "pi_fst_chr1_10kb.jpg"), plot = p_10kb, width = 7.5, height = 6, dpi = 300)

# --------------------------
# Chromosome I, 10 kb bins
p_10kb <- plot_pi_tajima_fst_combined(
  diversity_values_10kb,
  pi_ratios_10kb,
  fst_values_10kb,
  tajima_sample = "Treatment11_Pre",
  pi_sample = "Treatment11_Pre",
  fst_pair = "Treatment02_Pre:Treatment11_Pre",
  chrom = target_chrom,
  gene_region = tbb_iso1_region
) + plot_annotation(tag_levels = 'A')
ggsave(paste0(out_dir, "pi_tajima_fst_chr1_10kb.jpg"), plot = p_10kb, width = 7.5, height = 6, dpi = 300)
# --------------------------

# --------------------------
# Chromosome I, 100 kb bins
p_100kb <- plot_pi_tajima_fst_combined(
  diversity_values_100kb,
  pi_ratios_100kb,
  fst_values_100kb,
  tajima_sample = "Treatment11_Pre",
  pi_sample = "Treatment11_Pre",
  fst_pair = "Treatment02_Pre:Treatment11_Pre",
  chrom = target_chrom,
  gene_region = tbb_iso1_region
) + plot_annotation(tag_levels = 'A')
ggsave(paste0(out_dir, "pi_tajima_fst_chr1_100kb.jpg"), plot = p_100kb, width = 7.5, height = 6, dpi = 300)
# --------------------------
