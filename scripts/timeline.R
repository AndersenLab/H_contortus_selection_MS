library(tidyverse)
library(lubridate)
library(png)
library(patchwork)

# =============================================================================
# GLOBAL FONT SIZE SETTING
# =============================================================================
font_size <- 12

# Load clean data
timeline <- read_csv("~/Etta/H_contortus_selection/timeline_clean.csv")

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Format Dose_cat consistently for all data
timeline <- timeline %>%
  mutate(
    Dose_cat = case_when(
      Dose_cat == 0.25 ~ "0.25",
      Dose_cat == 0.5 ~ "0.50",
      Dose_cat == 0.75 ~ "0.75",
      Dose_cat == 5 ~ "5.00",
      TRUE ~ as.character(Dose_cat)
    ),
    Dose_cat = factor(Dose_cat, levels = c("0.25", "0.50", "0.75", "5.00")),
    Collection_Phase = factor(Collection_Phase, levels = c("Pre", "Post"))
  )

# Dataset for treatments
timeline_treatments <- timeline %>%
  filter(!is.na(Treatment)) %>%
  mutate(
    Date_label = format(Date, "%b %d, %Y"),
    label_nudge_x = case_when(
      Treatment == 10 ~ -5,
      Treatment == 11 ~ 5,
      TRUE ~ 0
    ),
    date_label_nudge_x = case_when(
      Treatment == 2 ~ -5,
      Treatment == 4 ~ 5,
      Treatment == 10 ~ -5,
      Treatment == 11 ~ 5,
      TRUE ~ 0
    )
  )

# Dataset for collections
timeline_collections <- timeline %>%
  filter(!is.na(Collection))

# Define colors
dose_colors <- c(
  "0.25" = viridisLite::viridis(10)[10],
  "0.50" = viridisLite::viridis(10)[7],
  "0.75" = viridisLite::viridis(10)[4],
  "5.00" = viridisLite::viridis(10)[2]
)

# Global min and max dates
date_min <- min(timeline$Date)
date_max <- max(timeline$Date)

# Y positions
y_treat <- 0.65
y_collect <- 0.25

# =============================================================================
# TIMELINE PLOT FUNCTION
# =============================================================================
make_timeline_plot <- function(color_collections_by_dose = FALSE) {
  
  p <- ggplot() +
    
    geom_segment(data = timeline_collections,
                 aes(x = Date, xend = Date, y = y_treat - 0.03, yend = y_collect + 0.03),
                 color = "grey70", linewidth = 0.3, linetype = "dotted") +
    
    geom_segment(aes(x = date_min, xend = date_max, y = y_treat, yend = y_treat),
                 color = "grey40", linewidth = 0.8) +
    
    geom_segment(aes(x = date_min, xend = date_max, y = y_collect, yend = y_collect),
                 color = "grey60", linewidth = 0.5) +
    
    geom_point(data = timeline_treatments, aes(x = Date, y = y_treat, color = Dose_cat),
               size = 4, shape = 16) +
    
    geom_text(data = timeline_treatments, 
              aes(x = Date + label_nudge_x, y = y_treat + 0.12, label = Treatment),
              size = font_size/.pt, color = "black") +
    
    annotate("text", x = date_min - 15, y = y_treat,
             label = "Treatment", size = font_size/.pt, color = "black", fontface = "bold", hjust = 1) +
    
    annotate("text", x = date_min - 15, y = y_collect,
             label = "Collection", size = font_size/.pt, color = "black", fontface = "bold", hjust = 1) +
    
    geom_text(data = timeline_treatments,
              aes(x = Date + date_label_nudge_x, y = y_treat + 0.28, label = Date_label),
              size = font_size/.pt, color = "black", angle = 45, hjust = 0, vjust = 0) +
    
    scale_color_manual(values = dose_colors, name = "Dose (mg/kg)") +
    scale_shape_manual(values = c("Pre" = 16, "Post" = 17), name = "Phase") +
    scale_x_date(expand = expansion(mult = c(0.02, 0.08))) +
    scale_y_continuous(limits = c(0.1, 1.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    theme_classic(base_size = font_size) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.spacing = unit(2, "pt"),
      legend.spacing.y = unit(2, "pt"),
      legend.title = element_text(size = font_size, color = "black"),
      legend.text = element_text(size = font_size, color = "black"),
      plot.margin = margin(t = 60, r = 50, b = 5, l = 80)
    )
  
  if(color_collections_by_dose) {
    p <- p + geom_point(data = timeline_collections, 
                        aes(x = Date, y = y_collect, shape = Collection_Phase, color = Dose_cat),
                        size = 2)
  } else {
    p <- p + geom_point(data = timeline_collections, 
                        aes(x = Date, y = y_collect, shape = Collection_Phase),
                        color = "grey50", size = 2)
  }
  
  return(p)
}

# =============================================================================
# CREATE AND SAVE STANDALONE TIMELINE PLOTS
# =============================================================================

p_timeline_grey <- make_timeline_plot(color_collections_by_dose = FALSE)
ggsave(
  filename = "~/Etta/H_contortus_selection/figures/timeline_collections_grey.png",
  plot = p_timeline_grey,
  width = 7.5,
  height = 3.3,
  dpi = 300
)

p_timeline_colored <- make_timeline_plot(color_collections_by_dose = TRUE)
ggsave(
  filename = "~/Etta/H_contortus_selection/figures/timeline_collections_colored.png",
  plot = p_timeline_colored,
  width = 7.5,
  height = 3.3,
  dpi = 300
)

# =============================================================================
# PADDOCK MAP PLOT
# =============================================================================
img <- readPNG("~/Etta/H_contortus_selection/Goat_fields.png")

img_height <- dim(img)[1]
img_width <- dim(img)[2]

p_paddock <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  annotate("text", x = 0.85, y = 0.75, 
           label = "Treated (6400 m²)", 
           size = font_size/.pt, color = "black", hjust = 1) +
  annotate("text", x = 0.85, y = 0.30, 
           label = "Control (2000 m²)", 
           size = font_size/.pt, color = "black", hjust = 1) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed(ratio = img_height/img_width, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 10))

# =============================================================================
# COMBINED FIGURE
# =============================================================================

p_timeline_grey_tight <- p_timeline_grey +
  theme(plot.margin = margin(t = 20, r = 50, b = 5, l = 80))

combined_plot_vertical <- p_paddock / p_timeline_grey_tight +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(size = font_size, face = "bold"))
  ) +
  plot_layout(heights = c(1.4, 0.8))

ggsave(
  filename = "~/Etta/H_contortus_selection/figures/figure1_combined.png",
  plot = combined_plot_vertical,
  width = 7.5,
  height = 7,
  dpi = 300
)
