# --------------------------
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(patchwork)
library(here)
source(here("scripts", "functions.R"))
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
# --------------------------

# --------------------------
# --- Comparison 1: both at 10kb ---
diversity_data_10kb <- load_diversity(here("data", "diversity", "10kb", "diversity.csv"))
diversity_values_10kb <- diversity_data_10kb$values
pi_ratios_10kb <- normalize_and_calculate_pi_ratio(diversity_values_10kb, baseline = "Treatment02_Pre")

fst_data_10kb <- load_fst(here("data", "fst", "10kb", "fst.csv"))
fst_values_10kb <- fst_data_10kb$values

# --- Comparison 2: both at 100kb ---
diversity_data_100kb <- load_diversity(here("data", "diversity", "100kb", "diversity.csv"))
diversity_values_100kb <- diversity_data_100kb$values
pi_ratios_100kb <- normalize_and_calculate_pi_ratio(diversity_values_100kb, baseline = "Treatment02_Pre")

fst_data_100kb <- load_fst(here("data", "fst", "100kb", "fst.csv"))
fst_values_100kb <- fst_data_100kb$values
# --------------------------

# --------------------------
# Reference region (isotype-1 beta-tubulin)
tbb_iso1_region <- c(7.027492, 7.031447)
target_chrom <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
# --------------------------

# --------------------------
# Define FST pairs: all comparisons to Treatment02_Pre baseline
fst_pairs_to_plot <- fst_values_100kb %>%
  filter(sample_1 == "Treatment02_Pre") %>%
  pull(sample_pair) %>%
  unique()
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
                                        smooth = FALSE,
                                        roman_numerals = FALSE,
                                        save_path = NULL,
                                        width = 7.5,
                                        height = 5,
                                        dpi = 300) {
  
  # Black color for points
  point_color <- "black"
  
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
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = NULL, y = "**log<sub>10</sub>(π ratio)**") +
    my_theme
  
  if (smooth) {
    p_pi <- p_pi + geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5)
  }
  
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
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = NULL, y = "**Tajima's D**") +
      my_theme
    
    if (smooth) {
      p_taj <- p_taj + geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5)
    }
    
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
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "**Position (Mb)**", y = "**F<sub>ST</sub>**") +
    my_theme
  
  if (smooth) {
    p_fst <- p_fst + geom_smooth(method = "loess", span = loess_span, se = FALSE, color = "black", linewidth = 0.5)
  }
  
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
out_dir <- here("figures")
# --------------------------

# --------------------------
# Chromosome I, 10 kb bins (pi + FST only)
p_10kb <- plot_pi_tajima_fst_combined(
  pi_ratios_10kb,
  fst_values_10kb,
  pi_sample = "Treatment11_Pre",
  fst_pair = "Treatment02_Pre:Treatment11_Pre",
  chrom = target_chrom,
  gene_region = tbb_iso1_region  
) + plot_annotation(tag_levels = 'A')

p_10kb
ggsave(file.path(out_dir, "pi_fst_chr1_10kb.jpg"), plot = p_10kb, width = 7.5, height = 6, dpi = 300)
# --------------------------

# --------------------------
# Define plot_fst_genome()
# Input: fst values, sample pair, optional gene region
# Output: Genome-wide FST plot faceted by chromosome

plot_fst_genome <- function(fst_values,
                            fst_pair,
                            gene_region = NULL,
                            gene_chrom = NULL,
                            compact = FALSE,
                            roman_numerals = FALSE,
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
                            height = 3,
                            dpi = 300) {
  
  # Black color for points
  point_color <- "black"
  
  # Chromosome labels
  if (roman_numerals) {
    chrom_labels <- c("I", "II", "III", "IV", "V", "X")
  } else {
    chrom_labels <- c("1", "2", "3", "4", "5", "X")
  }
  names(chrom_labels) <- chrom_order
  
  # X-axis breaks function
  x_breaks_fn <- if (compact) {
    function(x) setdiff(pretty(x), 0)
  } else {
    function(x) intersect(c(15, 30, 45), floor(min(x)):ceiling(max(x)))
  }
  
  # Filter and prepare data
  fst_data <- fst_values %>%
    dplyr::filter(sample_pair == fst_pair, !is.na(fst)) %>%
    dplyr::filter(chrom %in% chrom_order) %>%
    mutate(
      position_mb = ((start + end) / 2) / 1e6,
      chrom = factor(chrom, levels = chrom_order)
    )
  
  # Build plot
  p <- ggplot(fst_data, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.4, size = 0.6, color = point_color) +
    facet_grid(. ~ chrom, scales = "free_x", space = "free_x",
               labeller = labeller(chrom = chrom_labels)) +
    scale_x_continuous(breaks = x_breaks_fn, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "**Position (Mb)**", y = "**F<sub>ST</sub>**") +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = ggtext::element_markdown(face = "bold", size = 12),
      axis.title.y = ggtext::element_markdown(face = "bold", size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.1, "lines")
    )
  
  # Add gene marker if provided
  if (!is.null(gene_region) && !is.null(gene_chrom)) {
    gene_mid <- mean(gene_region)
    gene_data <- data.frame(
      chrom = factor(gene_chrom, levels = chrom_order),
      x = gene_mid
    )
    p <- p + geom_text(data = gene_data, aes(x = x, y = Inf, label = "*"),
                       vjust = 1.5, size = 8, color = "black", inherit.aes = FALSE)
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  p
}
# --------------------------

# --------------------------
# Genome-wide FST plot, 100kb windows
plot_fst_genome(
  fst_values = fst_values_100kb,
  fst_pair = "Treatment02_Pre:Treatment11_Pre",
  gene_region = tbb_iso1_region,
  gene_chrom = target_chrom,
  save_path = file.path(out_dir, "fst_genome_100kb.png")
)
# --------------------------

# --------------------------
# Define plot_fst_genome_multipair()
# Input: fst values, vector of sample pairs to plot
# Output: Genome-wide FST plot faceted by chromosome (columns) and sample pair (rows)

plot_fst_genome_multipair <- function(fst_values,
                                      fst_pairs,
                                      gene_region = NULL,
                                      gene_chrom = NULL,
                                      compact = FALSE,
                                      roman_numerals = FALSE,
                                      smooth = FALSE,
                                      loess_span = 0.1,
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
                                      height_per_pair = 1,
                                      dpi = 300) {
  
  # Set chromosome labels based on roman_numerals
  if (roman_numerals) {
    short_labels <- c("I", "II", "III", "IV", "V", "X")
  } else {
    short_labels <- c("1", "2", "3", "4", "5", "X")
  }
  
  # Filter and prepare data
  plot_data <- fst_values %>%
    dplyr::filter(sample_pair %in% fst_pairs, chrom %in% chrom_order, !is.na(fst)) %>%
    mutate(
      chrom = factor(chrom, levels = chrom_order),
      chrom_short = factor(short_labels[match(chrom, chrom_order)], levels = short_labels),
      position_mb = ((start + end) / 2) / 1e6
    )
  
  # Create human-readable labels from second sample in pair (sample_2)
  if (compact) {
    plot_data <- plot_data %>%
      mutate(sample_label = gsub("Treatment0?(\\d+)_(Pre|Post)", "\\2-Txt \\1", sample_2),
             sample_label = factor(sample_label, levels = unique(sample_label[order(match(sample_pair, fst_pairs))])))
  } else {
    plot_data <- plot_data %>%
      mutate(sample_label = gsub("Treatment0?(\\d+)_(Pre|Post)", "\\2-treatment \\1", sample_2),
             sample_label = factor(sample_label, levels = unique(sample_label[order(match(sample_pair, fst_pairs))])))
  }
  
  # Standard theme (matching plot_pi_ratio)
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
  
  # Set x-axis breaks (matching plot_pi_ratio)
  if (compact) {
    x_scale <- scale_x_continuous(breaks = function(x) setdiff(pretty(x), 0), expand = expansion(mult = c(0.02, 0.02)))
  } else {
    x_scale <- scale_x_continuous(breaks = c(15, 30, 45), expand = expansion(mult = c(0.02, 0.02)))
  }
  
  # Black color for points
  point_color <- "black"
  
  # Build plot
  p <- ggplot(plot_data, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.3, size = 0.6, color = point_color) +
    x_scale +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "Position (Mb)", y = "**F<sub>ST</sub>**") +
    my_theme +
    theme(axis.title.y = ggtext::element_markdown(size = 12))
  
  # Add loess smoothing if requested
  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = loess_span, se = FALSE, 
                         color = "black", linewidth = 0.5)
  }
  
  # Facet by chromosome (columns) and sample (rows)
  if (length(fst_pairs) == 1) {
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
  
  # Add gene marker if provided
  if (!is.null(gene_region) && !is.null(gene_chrom)) {
    gene_mid <- mean(gene_region)
    gene_chrom_short <- short_labels[match(gene_chrom, chrom_order)]
    
    # Create data for all sample labels at the gene chromosome
    gene_data <- plot_data %>%
      distinct(sample_label) %>%
      mutate(
        chrom_short = factor(gene_chrom_short, levels = short_labels),
        x = gene_mid
      )
    
    p <- p + geom_text(data = gene_data, aes(x = x, y = Inf, label = "*"),
                       vjust = 1.5, size = 8, color = "black", inherit.aes = FALSE)
  }
  
  if (!is.null(save_path)) {
    total_height <- height_per_pair * length(fst_pairs)
    ggsave(save_path, plot = p, width = width, height = total_height, dpi = dpi)
  }
  
  p
}
# --------------------------

# --------------------------
# Multi-pair FST plots (all timepoints vs Treatment02_Pre)

# Without smoothing (default) - Arabic numerals
plot_fst_genome_multipair(
  fst_values_100kb,
  fst_pairs = fst_pairs_to_plot,
  save_path = file.path(out_dir, "fst_genome_multipair_arabic.jpg")
)

# Without smoothing - Roman numerals
plot_fst_genome_multipair(
  fst_values_100kb,
  fst_pairs = fst_pairs_to_plot,
  roman_numerals = TRUE,
  save_path = file.path(out_dir, "fst_genome_multipair.jpg")
)

# With loess smoothing
plot_fst_genome_multipair(
  fst_values_100kb,
  fst_pairs = fst_pairs_to_plot,
  smooth = TRUE,
  loess_span = 0.1,
  save_path = file.path(out_dir, "fst_genome_multipair_smooth.jpg")
)
# --------------------------

