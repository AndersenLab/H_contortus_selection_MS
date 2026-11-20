library(tidyverse)
library(ggtext)
library(here)
source(here("scripts", "analysis", "functions.R"))

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

baseline <- "Treatment02_Pre"

plot_samples <- c("Treatment03_Pre", "Treatment04_Pre",
                  "Treatment08_Pre", "Treatment08_Post",
                  "Treatment09_Pre", "Treatment09_Post",
                  "Treatment10_Pre", "Treatment10_Post",
                  "Treatment11_Pre")

# -----------------------------------------------------------------------------
# Plotting function (unchanged from diversity_plots.R)
# -----------------------------------------------------------------------------

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
    x_scale <- scale_x_continuous(breaks = function(x) setdiff(pretty(x), 0),
                                  expand = expansion(mult = c(0.02, 0.02)))
  } else {
    x_scale <- scale_x_continuous(breaks = c(15, 30, 45),
                                  expand = expansion(mult = c(0.02, 0.02)))
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

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

diversity_values <- load_diversity(here("data", "diversity", "100kb", "diversity.csv"))$values

pi_ratios <- normalize_and_calculate_pi_ratio(diversity_values, baseline = baseline)

# -----------------------------------------------------------------------------
# Figure 2S
# -----------------------------------------------------------------------------

p_fig2s <- plot_pi_ratio(
  pi_ratios,
  samples = plot_samples,
  roman_numerals = FALSE
)

ggsave(
  here("figures", "figure2s_pi_ratio.jpg"),
  plot = p_fig2s,
  width = 7.5,
  height = 9,
  dpi = 300
)
