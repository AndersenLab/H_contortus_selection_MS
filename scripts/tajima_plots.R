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
# --------------------------
# 3. Plotting functions 
# --------------------------
# --------------------------

# --------------------------
# Define plot_tajima()
# Input: grenedalf dataframe with diversity metrics, which timepoints you want to plot
# Output: Genome-wide plot of Tajima D, faceted by sample, colored by chromosome

plot_tajima <- function(data,
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
                        show_label = FALSE,
                        smooth = FALSE,
                        loess_span = 0.1,
                        save_path = NULL,
                        width = 10,
                        height = 3,
                        dpi = 300) {
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(tajimas_d)) %>%
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
  
  p <- ggplot(plot_data, aes(x = position_mb, y = tajimas_d, color = chrom_short)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_viridis_d() +
    scale_x_continuous(breaks = axis_labels$mid, labels = axis_labels$chrom_short) +
    labs(x = "Chromosome", y = "Tajima D", color = "Chromosome") +
    guides(color = "none") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  # Add facet with or without labels
  if (show_label) {
    p <- p + 
      facet_wrap(~sample, ncol = 1, strip.position = "left", labeller = as_labeller(relabel_sample)) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)
      )
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
# Define plot_tajima_by_value()
# Input: diversity dataframe, one or more samples to plot
# Output: Genome-wide plot colored by Tajima D value instead of chromosome

plot_tajima_by_value <- function(data,
                                 samples,
                                 chrom_order,
                                 save_path = NULL,
                                 width = 10,
                                 height = 3,
                                 dpi = 300) {
  
  plot_data <- data %>%
    dplyr::filter(sample %in% samples, chrom %in% chrom_order, !is.na(tajimas_d)) %>%
    mutate(position_mb = (start + end) / 2 / 1e6)
  
  p <- ggplot(plot_data, aes(x = position_mb, y = tajimas_d, color = tajimas_d)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_viridis_c() +
    facet_wrap(~sample, ncol = 1, strip.position = "left", labeller = as_labeller(relabel_sample)) +
    labs(x = "Position (Mb)", y = "Tajima D", color = "Tajima D") +
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

# Get all sample names for plotting
all_samples <- unique(diversity_values$sample)
# --------------------------

# --------------------------
# --------------------------
# 5. Generate plots 
# --------------------------
# --------------------------

# Define output directory
out_dir <- "/vast/eande106/projects/Etta/H_contortus_selection/figures/"

# --------------------------
# Plot Tajima D for single sample
plot_tajima(
  diversity_values,
  samples = "Treatment11_Pre"
)
ggsave(paste0(out_dir, "tajima_Treatment11_Pre.jpg"), plot = last_plot(), width = 3.75, height = 2.5, dpi = 300)
# --------------------------

# --------------------------
# Plot Tajima D for all samples
plot_tajima(
  diversity_values,
  samples = c("Treatment03_Pre", "Treatment04_Pre", "Treatment08_Pre",
              "Treatment08_Post", "Treatment09_Pre", "Treatment09_Post",
              "Treatment10_Pre", "Treatment10_Post", "Treatment11_Pre"),
  show_label = TRUE
)
ggsave(paste0(out_dir, "tajima_all_samples.jpg"), plot = last_plot(), width = 7.5, height = 9, dpi = 300)
# --------------------------

# --------------------------
# Plot Tajima D colored by value (single chromosome)
plot_tajima_by_value(
  diversity_values,
  samples = all_samples,
  chrom_order = "hcontortus_chr1_Celeg_TT_arrow_pilon"
)
ggsave(paste0(out_dir, "tajima_by_value_chr1.jpg"), plot = last_plot(), width = 3.75, height = 8, dpi = 300)
# --------------------------