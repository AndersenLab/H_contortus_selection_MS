# --------------------------
# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtext)
# --------------------------

# --------------------------
# Set working directory 
working_directory <- "~/Etta/H_contortus_selection_paper/"
setwd(working_directory)
# --------------------------

# --------------------------
# Define load_fst()
# Input: Raw FST output from grenedalf (fst.csv)
# Output: A list of three dataframes
# 1. windows: The number of sites that were included or excluded by window 
# 2. qc: The number of sites that were included or excluded by window and by sample pair
# 3. values: The numeric output (just FST value) by window and by sample pair

load_fst <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE) # Read raw data, suppress warning messages

  # Select columns for windows dataframe
  windows <- raw %>%
    dplyr::select(chrom, start, end, starts_with("total."))

  # Each comparison is a column (Example: Treatment09_Pre:Treatment09_Post.passed)
  # Convert to long format where each comparison gets its own row
  long <- raw %>%
    dplyr::select(chrom, start, end, matches("^[^.]+:[^.]+\\.(missing|numeric|passed|fst)$")) %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = c("sample_pair", "metric"),
      names_pattern = "^(.+)\\.(missing|numeric|passed|fst)$",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value)) #Convert NaN to NA 
  
  # Make columns for individual timepoint labels 
  # (Example: For Treatment09_Pre:Treatment09_Post, Treatment09_Pre is sample 1, Treatment09_Post is sample 2)
  wide <- long %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    separate(
      sample_pair,
      into = c("sample_1", "sample_2"),
      sep = ":",
      remove = FALSE
    )
  
  # Make qc dataframe
  qc <- wide %>%
    dplyr::select(chrom, start, end, sample_pair, sample_1, sample_2, missing, numeric, passed)
  
  # Make fst values dataframe
  values <- wide %>%
    dplyr::select(chrom, start, end, sample_pair, sample_1, sample_2, fst)
  
  # Store as a list
  list(
    windows = windows,
    qc = qc,
    values = values
  )
}
# --------------------------

# --------------------------
# Load and process raw data with load_fst()
fst_data <- load_fst("10X/fst.csv")

fst_windows <- fst_data$windows
fst_qc <- fst_data$qc
fst_values <- fst_data$values
# --------------------------

# --------------------------
# Prepare data for plotting

# Make a list of chromosome names
chrom_ids <- c(
  "hcontortus_chr1_Celeg_TT_arrow_pilon",
  "hcontortus_chr2_Celeg_TT_arrow_pilon",
  "hcontortus_chr3_Celeg_TT_arrow_pilon",
  "hcontortus_chr4_Celeg_TT_arrow_pilon",
  "hcontortus_chr5_Celeg_TT_arrow_pilon",
  "hcontortus_chrX_Celeg_TT_arrow_pilon"
) 

# Look for integers or X after "chr"
chrom_num <- str_extract(chrom_ids, "(?<=chr)[0-9X]+") 

# Convert integers to roman numerals
roman_labs <- ifelse(chrom_num == "X", "X", as.character(as.roman(as.integer(chrom_num)))) 

# Make list of nice labels for plotting
chrom_labels <- setNames(
  paste0("*H. contortus* chromosome ", roman_labs),
  chrom_ids
)
# --------------------------

# --------------------------
# Define plot_fst()
# Input: data (fst values from load_fst()), and pair (timepoint comparison)
# Output: Manhattan-style-plot for all chromosomes

# Define function
plot_fst <- function(data,
                     pair,
                     chroms = NULL,
                     exclude_chroms = "mitochondrion",
                     vline_pos = NULL,
                     gene_region = NULL,
                     smooth = TRUE,
                     loess_span = 0.1,
                     labels = chrom_labels,
                     save_path = NULL,
                     width = 8,
                     height = 5,
                     dpi = 300) {
  plot_data <- data %>%
    dplyr::filter(sample_pair == pair, !is.na(fst)) # Filter to selected pair, drop NAs
  
# If chroms is supplied with a value, just plot that value. Otherwise do everything but mitochondria
  if (!is.null(chroms)) {
    plot_data <- plot_data %>% dplyr::filter(chrom %in% chroms)
  } else if (!is.null(exclude_chroms)) {
    plot_data <- plot_data %>% dplyr::filter(!chrom %in% exclude_chroms)
  }
  
# Calculate midpoint for each window
  
  plot_data <- plot_data %>%
    mutate(position_mb = ((start + end) / 2) / 1e6)

# If the chromosome is not in the nice labels, use the raw name
  facet_labeller <- as_labeller(function(x) ifelse(x %in% names(labels), labels[x], x))
  
# Extract timepoint pair names for the title
  pair_label <- paste(unique(plot_data$sample_1), "vs", unique(plot_data$sample_2))
  
# Make the plot
  p <- ggplot(plot_data, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.3, size = 0.8) +
    facet_wrap(~chrom, scales = "free_x", labeller = facet_labeller) +
    labs(
      x = "Position (Mb)", y = expression(F[ST]),
      title = pair_label
    ) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    theme_bw() +
    theme(strip.text = element_markdown())

# Annotate the plot with region of interest
  if (!is.null(gene_region)) {
    p <- p + annotate("rect",
      xmin = gene_region[1], xmax = gene_region[2],
      ymin = -Inf, ymax = Inf,
      fill = "red", alpha = 0.15
    )

    gene_midpoint <- mean(gene_region)
    p <- p + geom_vline(
      xintercept = gene_midpoint, linetype = "dashed",
      color = "red", linewidth = 0.6
    )
  }
  
# Smooth with local regression
  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = loess_span, se = FALSE)
  }
  
# Draw dashed line at isotype-1 or isotype-2
  
  if (!is.null(vline_pos)) {
    p <- p + geom_vline(xintercept = vline_pos, linetype = "dashed", color = "red")
  }
  
# Save plot

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
  }

  p
}
# --------------------------

# --------------------------
# Save locations of isotype-1 and isotype-2 beta-tubulins
tbb_iso1_region <- c(7.027492, 7.031447)
tbb_iso2_region <- c(13.433191, 13.438615)
# --------------------------

# --------------------------
# Make plots by chromosome

fst_data <- load_fst("10X/fst.csv")
fst_values <- fst_data$values

plot_fst(
  fst_values,
  pair = "Treatment02_Pre:Treatment11_Pre",
  chroms = "hcontortus_chr1_Celeg_TT_arrow_pilon",
  gene_region = tbb_iso1_region,
  # save_path = "fst_chr1_earliest_vs_latest.png"   # omit save_path to just view in-session
)

ggsave("chrX_fst.png", plot = last_plot(), width = 8, height = 5, dpi = 300, scale = 0.7)
# --------------------------

# --------------------------
# Define plot_fst_horizontal()
# Makes one plot for the whole genome 

plot_fst_horizontal <- function(data,
                                pair,
                                chroms = NULL,
                                exclude_chroms = "mitochondrion",
                                vline_pos = NULL,
                                gene_region = NULL,
                                gene_chrom = NULL,
                                marker_style = "asterisk",
                                smooth = TRUE,
                                loess_span = 0.1,
                                labels = chrom_labels,
                                save_path = NULL,
                                width = 7.5,
                                height = 2.5,
                                dpi = 300) {
  plot_data <- data %>%
    dplyr::filter(sample_pair == pair, !is.na(fst))

  if (!is.null(chroms)) {
    plot_data <- plot_data %>% dplyr::filter(chrom %in% chroms)
  } else if (!is.null(exclude_chroms)) {
    plot_data <- plot_data %>% dplyr::filter(!chrom %in% exclude_chroms)
  }

  plot_data <- plot_data %>%
    mutate(position_mb = ((start + end) / 2) / 1e6)

  # Define chromosome order
  chrom_order <- c(
    "hcontortus_chr1_Celeg_TT_arrow_pilon",
    "hcontortus_chr2_Celeg_TT_arrow_pilon",
    "hcontortus_chr3_Celeg_TT_arrow_pilon",
    "hcontortus_chr4_Celeg_TT_arrow_pilon",
    "hcontortus_chr5_Celeg_TT_arrow_pilon",
    "hcontortus_chrX_Celeg_TT_arrow_pilon"
  )

  plot_data <- plot_data %>%
    dplyr::filter(chrom %in% chrom_order) %>%
    mutate(chrom = factor(chrom, levels = chrom_order))

  # Drop full chromosome names, just label with Roman numerals
  short_labels <- setNames(
    c("I", "II", "III", "IV", "V", "X"),
    chrom_order
  )

  facet_labeller <- as_labeller(function(x) short_labels[x])

  p <- ggplot(plot_data, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.3, size = 0.4) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x", labeller = facet_labeller) +
    labs(x = "Position (Mb)", y = expression(F[ST])) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 9),
      panel.spacing = unit(0.3, "lines")
    )

  # Gene region highlighting
  if (!is.null(gene_region) && !is.null(gene_chrom)) {
    gene_df <- data.frame(
      chrom = factor(gene_chrom, levels = chrom_order),
      x = mean(gene_region),
      y = 0.6
    )

    if (marker_style == "asterisk") {
      p <- p + geom_text(
        data = gene_df, aes(x = x, y = y),
        label = "*", size = 6, color = "red"
      )
    } else if (marker_style == "line") {
      p <- p + geom_vline(
        data = gene_df, aes(xintercept = x),
        linetype = "dashed", color = "red", linewidth = 0.5
      )
    }
  }

  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = loess_span, se = FALSE, linewidth = 0.7)
  }

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
  }

  p
}
# --------------------------

# --------------------------

# Make plots
# With asterisk
plot_fst_horizontal(
  fst_values,
  pair = "Treatment02_Pre:Treatment11_Pre",
  gene_region = tbb_iso1_region,
  gene_chrom = "hcontortus_chr1_Celeg_TT_arrow_pilon",
  marker_style = "asterisk",
  save_path = "fst_asterisk.png"
)

# With dashed line
plot_fst_horizontal(
  fst_values,
  pair = "Treatment02_Pre:Treatment11_Pre",
  gene_region = tbb_iso1_region,
  gene_chrom = "hcontortus_chr1_Celeg_TT_arrow_pilon",
  marker_style = "line",
  save_path = "fst_line.png"
)
# --------------------------
