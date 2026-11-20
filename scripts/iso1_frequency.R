library(tidyverse)
library(ggtext)
library(patchwork)
library(here)

# ---- Shared constants -------------------------------------------------

timepoint_levels <- c("1","2","3","4","5","6","7","8","9","10","11")
phase_levels     <- c("Pre", "Post")

# ---- Data processing: sync file (position 7029790) --------------------

sample_names <- c(
  "Treatment01_Pre", "Treatment01_Post",
  "Treatment02_Pre", "Treatment02_Post",
  "Treatment03_Pre", "Treatment03_Post",
  "Treatment04_Pre", "Treatment04_Post",
  "Treatment05_Pre", "Treatment05_Post",
  "Treatment06_Pre", "Treatment06_Post",
  "Treatment07_Pre", "Treatment07_Post",
  "Treatment08_Pre", "Treatment08_Post",
  "Treatment09_Pre", "Treatment09_Post",
  "Treatment10_Pre", "Treatment10_Post",
  "Treatment11_Pre", "Treatment11_Post"
)

sync_file <- here("data", "iso-1", "tbb_iso1_resistance_sites.sync")

sync_line <- read.table(
  sync_file,
  header = FALSE, sep = "\t", stringsAsFactors = FALSE
) %>%
  filter(V2 == 7029790)

counts <- sync_line %>%
  select(4:25) %>%
  setNames(sample_names) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Counts") %>%
  separate(Counts, into = c("A", "T", "C", "G", "N", "del"),
           sep = ":", convert = TRUE) %>%
  mutate(
    F200Y_count = A,
    REF_count   = T,
    Total       = F200Y_count + REF_count,
    F200Y_frequency = ifelse(Total > 0, 100 * F200Y_count / Total, 0),
    REF_frequency   = ifelse(Total > 0, 100 * REF_count   / Total, 0)
  )

f200y_sync <- counts %>%
  mutate(
    Treatment    = str_extract(Sample, "Treatment[0-9]+"),
    Phase        = str_extract(Sample, "Pre|Post"),
    TreatmentNum = as.numeric(str_remove(Treatment, "Treatment")),
    Timepoint    = as.character(TreatmentNum)
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    Phase     = factor(Phase, levels = phase_levels)
  ) %>%
  select(Timepoint, TreatmentNum, Phase, F200Y_frequency, REF_frequency, Total)

# ---- Bar plot data prep -----------------------------------------------

f200y_plot_data <- f200y_sync %>%
  mutate(
    x_combined = paste(Timepoint, Phase, sep = "_"),
    x_combined = factor(x_combined, levels = paste(rep(timepoint_levels, each = 2), 
                                                   rep(phase_levels, times = 11), 
                                                   sep = "_")),
    x_numeric = as.numeric(x_combined)
  )

depth_labels <- f200y_plot_data %>%
  select(x_combined, x_numeric, Phase, Timepoint, Total) %>%
  distinct() %>%
  mutate(Label = paste0("n=", Total))

timepoint_boxes <- tibble(
  Timepoint = factor(timepoint_levels, levels = timepoint_levels),
  xmin = seq(0.5, by = 2, length.out = 11) + 0.1,
  xmax = seq(2.5, by = 2, length.out = 11) - 0.1
)

# ---- Bar plot function ------------------------------------------------

make_f200y_barplot <- function(df, depth_labels, timepoint_boxes) {
  
  df_long <- df %>%
    pivot_longer(cols = c(F200Y_frequency, REF_frequency),
                 names_to = "Category", values_to = "Percent") %>%
    mutate(Category = factor(Category,
                             levels = c("REF_frequency", "F200Y_frequency"),
                             labels = c("Susceptible (REF)", "Resistant (F200Y)")))
  
  variant_colors <- c(
    "Susceptible (REF)" = "blue",
    "Resistant (F200Y)" = "red"
  )
  
  x_labels <- rep(c("Pre", "Post"), 11)
  
  p <- ggplot(df_long, aes(x = x_numeric, y = Percent, fill = Category)) +
    geom_rect(data = timepoint_boxes,
              aes(xmin = xmin, xmax = xmax, ymin = -15, ymax = -5),
              inherit.aes = FALSE,
              fill = "grey85", color = NA) +
    geom_text(data = timepoint_boxes,
              aes(x = (xmin + xmax) / 2, y = -10, label = Timepoint),
              inherit.aes = FALSE,
              size = 12 / .pt) +
    geom_col(position = position_stack(reverse = FALSE), width = 0.7) +
    geom_text(data = depth_labels,
              aes(x = x_numeric, y = 106, label = Label),
              inherit.aes = FALSE,
              size = 12 / .pt,
              angle = 60,
              hjust = 0,
              vjust = 0.5) +
    scale_fill_manual(
      values = variant_colors,
      breaks = c("Resistant (F200Y)", "Susceptible (REF)")
    ) +
    scale_x_continuous(
      breaks = 1:22,
      labels = x_labels,
      expand = expansion(add = 0.5)
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 25),
      expand = expansion(mult = c(0, 0.15))
    ) +
    coord_cartesian(ylim = c(-15, 115), clip = "off") +
    labs(
      x = "**Treatment**",
      y = "**F200Y variant frequency (%)**",
      fill = "Genotype",
      tag = "A"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title.x = element_markdown(size = 12),
      axis.title.y = element_markdown(size = 12),
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt"),
      plot.tag = element_text(size = 14, face = "bold")
    )
  
  p
}

# ---- Scatter plot function --------------------------------------------

make_f200y_scatter <- function(df, min_coverage = 3) {
  
  df_filtered <- df %>%
    filter(Total >= min_coverage)
  
  p <- ggplot(df_filtered, aes(x = TreatmentNum, y = F200Y_frequency)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, fullrange = TRUE) +
    geom_point(aes(shape = Phase, size = Total), color = "red", alpha = 0.8) +
    scale_shape_manual(
      values = c("Pre" = 16, "Post" = 17),
      labels = c("Pre" = "Pre-treatment", "Post" = "Post-treatment"),
      name = "Phase"
    ) +
    scale_size_continuous(name = "Coverage", range = c(2, 6)) +
    scale_x_continuous(
      breaks = 1:11,
      labels = as.character(1:11),
      limits = c(0.5, 11.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 25)
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      x = "**Treatment**",
      y = "**F200Y variant frequency (%)**",
      tag = "B"
    ) +
    guides(
      shape = guide_legend(order = 1),
      size = guide_legend(order = 2)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title.x = element_markdown(size = 12),
      axis.title.y = element_markdown(size = 12),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 10, r = 5, b = 5, l = 5, unit = "pt"),
      plot.tag = element_text(size = 14, face = "bold")
    )
  
  p
}

# ---- Generate plots ---------------------------------------------------

plot_bar <- make_f200y_barplot(f200y_plot_data, depth_labels, timepoint_boxes)
plot_scatter <- make_f200y_scatter(f200y_sync, min_coverage = 5)

# ---- Combine with patchwork -------------------------------------------

combined_plot <- (plot_bar / plot_scatter) +
  plot_layout(axes = "collect_x")

combined_plot

# ---- Save -------------------------------------------------------------

dir.create(here("figures"), showWarnings = FALSE)

ggsave(
  here("figures", "figure4_F200Y_combined.jpg"),
  combined_plot,
  width = 7.5,
  height = 8,
  dpi = 300
)