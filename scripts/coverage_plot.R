library(ggplot2)
library(ggtext)
library(cowplot)

coverage <- read.delim("~/Etta/H_contortus_selection/coverage/coverage_with_names.tsv")

# Make names human readable
coverage$sample_name <- gsub("Treatment0?(\\d+)_(Pre|Post)", "\\2-treatment \\1", coverage$sample_name)

# Make timepoints factors so they are ordered correctly 
correct_order <- c(rbind(
  paste0("Pre-treatment ", 1:11),
  paste0("Post-treatment ", 1:11)
))

coverage$sample_name <- factor(coverage$sample_name, levels = rev(correct_order))

# Plot
# Create a shared theme
my_theme <- theme_bw(base_size = 12) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank())

# Horizontal version
p_horizontal <- ggplot(coverage, aes(x = coverage, y = sample_name)) + 
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 10, linetype = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlab("Mean coverage") +
  ylab("Timepoint") +
  my_theme

ggsave("~/Etta/H_contortus_selection/figures/coverage_horizontal.png", p_horizontal, width = 7.5, height = 4.5, dpi = 300)

# Flipped version
p_vertical <- ggplot(coverage, aes(x = sample_name, y = coverage)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 10, linetype = 3) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlab("Timepoint") +
  ylab("Mean coverage") +
  my_theme +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(5, 10, 5, 25))

ggsave("~/Etta/H_contortus_selection/figures/coverage_vertical.png", p_vertical, width = 7.5, height = 3.75, dpi = 300)
