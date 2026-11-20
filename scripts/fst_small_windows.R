#!/usr/bin/env Rscript
# =========================================================
# plot_fst.R
#
# Plots FST from a single grenedalf fst.csv, at two zoom levels:
#   1. Genome-wide, faceted by chromosome (Mb scale)
#   2. Zoomed in on a candidate gene region (bp scale), with known
#      resistance-variant sites marked
#
# Both plots pull from the SAME fst.csv -- so if you ran grenedalf
# genome-wide with small windows, you don't need a second run to get
# the zoomed view; this script just filters/re-scales the same data.
#
# Interactive use: run this in an R console / RStudio, section by section.
# Load the functions (top of the file) first, then work through the
# "PLOTTING" section at the bottom -- set fst_csv_path and pair, then
# run each plot_fst_zoom()/plot_fst() call as needed.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggtext)   # for italicized species name in facet strip text
  library(data.table)  # for fast, memory-efficient pre-filtered reads
})

# ---------------------------------------------------------
# load_fst(): reads a grenedalf fst.csv and returns three tidy tables
#
# Returns a list with:
#   windows - one row per chrom/start/end, the window-level total.* QC columns
#   qc      - one row per chrom/start/end/sample_pair, the per-pair
#             .missing/.numeric/.passed columns
#   values  - one row per chrom/start/end/sample_pair, just the .fst value
#
# sample_pair is kept as-is (e.g. "Treatment02_Pre:Treatment03_Pre") and also
# split into sample_1 / sample_2 in both qc and values.
#
# chroms: optional character vector of chromosome names to keep, e.g.
#   c("hcontortus_chr1_Celeg_TT_arrow_pilon", "hcontortus_chr2_Celeg_TT_arrow_pilon")
# When provided, filtering happens via awk BEFORE the file is read into R,
# so rows from other chromosomes are never loaded into memory at all.
#
# pair: optional single sample_pair string (e.g. "Treatment02_Pre:Treatment3_Post").
# When provided, only that pair's 4 columns (.missing/.numeric/.passed/.fst)
# are read, alongside chrom/start/end/total.* -- every other pair's columns
# are dropped before the file is read into R. Combined with `chroms`, this
# is the key trick for handling large genome-wide, fine-window fst.csv
# files with many sample pairs: you only ever load the rows and columns
# you actually need for one plot.
# ---------------------------------------------------------
load_fst <- function(path, chroms = NULL, pair = NULL) {
  
  path <- path.expand(path)  # resolve ~ before handing the path to awk/fread's cmd=
  
  needs_filter <- !is.null(chroms) || !is.null(pair)
  
  if (needs_filter) {
    header <- names(fread(path, nrows = 0))
    
    if (!is.null(pair)) {
      keep_cols <- c("chrom", "start", "end",
                     header[startsWith(header, "total.")],
                     paste0(pair, c(".missing", ".numeric", ".passed", ".fst")))
      missing_cols <- setdiff(paste0(pair, c(".missing", ".numeric", ".passed", ".fst")), header)
      if (length(missing_cols) > 0) {
        stop(sprintf(
          "Pair '%s' not found in %s. Check the exact pair name (see the csv header).",
          pair, path
        ))
      }
    } else {
      keep_cols <- header
    }
    
    col_idx <- match(keep_cols, header)
    field_list <- paste(paste0("$", col_idx), collapse = ",")
    
    row_filter <- "1"  # awk true, i.e. keep all rows by default
    if (!is.null(chroms)) {
      row_filter <- paste0('(', paste0('$1=="', chroms, '"', collapse = " || "), ')')
    }
    
    awk_cmd <- sprintf(
      'awk -F, -v OFS=, \'NR==1 { print %s; next } %s { print %s }\' %s',
      field_list, row_filter, field_list, shQuote(path)
    )
    
    message(
      "Pre-filtering", 
      if (!is.null(chroms)) paste0(" to chromosome(s): ", paste(chroms, collapse = ", ")) else "",
      if (!is.null(pair)) paste0(" | pair: ", pair) else ""
    )
    raw <- fread(cmd = awk_cmd, showProgress = FALSE) %>% as_tibble()
  } else {
    raw <- read_csv(path, show_col_types = FALSE)
  }
  
  windows <- raw %>%
    select(chrom, start, end, starts_with("total."))
  
  long <- raw %>%
    select(chrom, start, end, matches("^[^.]+:[^.]+\\.(missing|numeric|passed|fst)$")) %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = c("sample_pair", "metric"),
      names_pattern = "^(.+)\\.(missing|numeric|passed|fst)$",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(is.nan(value), NA_real_, value))
  
  wide <- long %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    separate(
      sample_pair,
      into = c("sample_1", "sample_2"),
      sep = ":",
      remove = FALSE
    )
  
  qc <- wide %>%
    select(chrom, start, end, sample_pair, sample_1, sample_2, missing, numeric, passed)
  
  values <- wide %>%
    select(chrom, start, end, sample_pair, sample_1, sample_2, fst)
  
  list(
    windows = windows,
    qc = qc,
    values = values
  )
}

# ---------------------------------------------------------
# chrom_labels: maps raw chrom IDs to display labels (Markdown italics
# via ggtext::element_markdown()). Add more entries here as needed.
# ---------------------------------------------------------
chrom_ids <- c("hcontortus_chr1_Celeg_TT_arrow_pilon",
               "hcontortus_chr2_Celeg_TT_arrow_pilon",
               "hcontortus_chr3_Celeg_TT_arrow_pilon",
               "hcontortus_chr4_Celeg_TT_arrow_pilon",
               "hcontortus_chr5_Celeg_TT_arrow_pilon",
               "hcontortus_chrX_Celeg_TT_arrow_pilon")

chrom_num <- str_extract(chrom_ids, "(?<=chr)[0-9X]+")
roman_labs <- ifelse(chrom_num == "X", "X", as.character(as.roman(as.integer(chrom_num))))

chrom_labels <- setNames(
  paste0("*H. contortus* chromosome ", roman_labs),
  chrom_ids
)

# ---------------------------------------------------------
# plot_fst(): genome-position Fst plot for one sample_pair, faceted by
# chromosome, x-axis in Mb. For a full genome-wide overview.
# ---------------------------------------------------------
plot_fst <- function(data,
                     pair,
                     chroms = NULL,
                     exclude_chroms = "mitochondrion",
                     vline_pos = NULL,
                     gene_region = NULL,      # in Mb
                     smooth = TRUE,
                     loess_span = 0.1,
                     labels = chrom_labels,
                     save_path = NULL,
                     width = 8,
                     height = 5,
                     dpi = 300) {
  
  plot_data <- data %>%
    filter(sample_pair == pair, !is.na(fst))
  
  if (nrow(plot_data) == 0) {
    stop(sprintf("No FST data found for sample_pair '%s'. Check the exact pair name in your fst.csv header.", pair))
  }
  
  if (!is.null(chroms)) {
    plot_data <- plot_data %>% filter(chrom %in% chroms)
  } else if (!is.null(exclude_chroms)) {
    plot_data <- plot_data %>% filter(!chrom %in% exclude_chroms)
  }
  
  plot_data <- plot_data %>%
    mutate(position_mb = ((start + end) / 2) / 1e6)
  
  facet_labeller <- as_labeller(function(x) ifelse(x %in% names(labels), labels[x], x))
  
  pair_label <- paste(unique(plot_data$sample_1), "vs", unique(plot_data$sample_2))
  
  p <- ggplot(plot_data, aes(x = position_mb, y = fst)) +
    geom_point(alpha = 0.3, size = 0.8) +
    facet_wrap(~ chrom, scales = "free_x", labeller = facet_labeller) +
    labs(x = "Position (Mb)", y = expression(F[ST]),
         title = pair_label) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    theme_bw() +
    theme(strip.text = element_markdown())
  
  if (!is.null(gene_region)) {
    p <- p + annotate("rect",
                      xmin = gene_region[1], xmax = gene_region[2],
                      ymin = -Inf, ymax = Inf,
                      fill = "red", alpha = 0.15)
    
    gene_midpoint <- mean(gene_region)
    p <- p + geom_vline(xintercept = gene_midpoint, linetype = "dashed",
                        color = "red", linewidth = 0.6)
  }
  
  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = loess_span, se = FALSE)
  }
  
  if (!is.null(vline_pos)) {
    p <- p + geom_vline(xintercept = vline_pos, linetype = "dashed", color = "red")
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
    message("Saved plot to: ", save_path)
  }
  
  p
}

# ---------------------------------------------------------
# plot_fst_zoom(): fine-scale Fst plot for a chromosome or region.
# X-axis is in bp (not Mb). No faceting -- filters down to one
# chromosome from a (possibly multi-chromosome) fst dataset.
#
# `filter_region` controls what data is DISPLAYED (i.e. what's kept
# after filtering) -- pass c(start_bp, end_bp) to zoom into just that
# span, or leave NULL to show the whole chromosome.
#
# `gene_region` controls what gets HIGHLIGHTED (the shaded rectangle) --
# independent of filter_region, so you can show a gene's location on a
# whole-chromosome plot without restricting the data to just the gene.
# ---------------------------------------------------------
plot_fst_zoom <- function(data,
                          pair,
                          chrom = NULL,
                          filter_region = NULL,       # c(start_bp, end_bp) -- restricts DISPLAYED data
                          flank = 0,                   # extra bp shown on each side of filter_region
                          gene_region = NULL,          # c(start_bp, end_bp) -- shaded rectangle only, doesn't filter
                          candidate_sites = NULL,      # named vector, e.g. c(F200Y = 7029790)
                          smooth = TRUE,
                          loess_span = 0.3,
                          title = NULL,
                          save_path = NULL,
                          width = 8,
                          height = 5,
                          dpi = 300) {
  
  plot_data <- data %>%
    filter(sample_pair == pair, !is.na(fst))
  
  if (!is.null(chrom)) {
    plot_data <- plot_data %>% filter(chrom == !!chrom)
  }
  
  if (!is.null(filter_region)) {
    lo <- filter_region[1] - flank
    hi <- filter_region[2] + flank
    plot_data <- plot_data %>%
      filter((start + end) / 2 >= lo, (start + end) / 2 <= hi)
  }
  
  if (nrow(plot_data) == 0) {
    stop(sprintf(
      "No FST data found for sample_pair '%s' in the specified region. Check the pair name, chrom, and filter_region.",
      pair
    ))
  }
  
  plot_data <- plot_data %>%
    mutate(position_bp = (start + end) / 2)
  
  pair_label <- paste(unique(plot_data$sample_1), "vs", unique(plot_data$sample_2))
  
  p <- ggplot(plot_data, aes(x = position_bp, y = fst)) +
    geom_point(alpha = 0.5, size = 1) +
    labs(x = "Position (bp)", y = expression(F[ST]),
         title = if (!is.null(title)) title else pair_label) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    theme_bw()
  
  if (!is.null(gene_region)) {
    p <- p + annotate("rect",
                      xmin = gene_region[1], xmax = gene_region[2],
                      ymin = -Inf, ymax = Inf,
                      fill = "red", alpha = 0.08)
  }
  
  if (!is.null(candidate_sites)) {
    site_df <- data.frame(pos = candidate_sites, label = names(candidate_sites))
    p <- p +
      geom_vline(data = site_df, aes(xintercept = pos),
                 linetype = "dashed", color = "darkred", linewidth = 0.5) +
      geom_text(data = site_df, aes(x = pos, y = 0.95, label = label),
                angle = 90, vjust = -0.4, hjust = 1, size = 3, color = "darkred")
  }
  
  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = loess_span, se = FALSE)
  }
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
    message("Saved plot to: ", save_path)
  }
  
  p
}

# ---------------------------------------------------------
# Known candidate loci and resistance-variant sites, in bp, on this genome
# assembly (arrow_pilon). Source: Doyle et al. chromosome-scale H. contortus
# assembly paper (isotype-1/isotype-2 gene coordinates), and this session's
# codon-position calculations from the WBPS19 GFF (resistance-variant sites).
# ---------------------------------------------------------
tbb_iso1_chrom  <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
tbb_iso1_region <- c(7027492, 7031447)
tbb_iso1_sites  <- c(
  F167Y = 7029569,
  E198A = 7029784,   # 2nd base of codon 198 -- also covers E198V etc.
  F200Y = 7029790
)

tbb_iso2_chrom  <- "hcontortus_chr2_Celeg_TT_arrow_pilon"
tbb_iso2_region <- c(13433191, 13438615)
tbb_iso2_sites  <- c(
  E198A = 13435823    # 2nd base of codon 198
)

# =========================================================
# PLOTTING -- interactive
#
# Set fst_csv_path and pair to match your data, then run the sections
# below one at a time. Re-run load_fst() whenever you point at a new
# fst.csv; the three plot calls all reuse fst_values.
# =========================================================

fst_csv_path <- "~/vast-eande106/projects/Etta/H_contortus_selection_paper/scripts/small_windows/fst.csv"                              # <- set this
pair         <- "Treatment02_Pre:Treatment10_Post"       # <- set this

# ---------------------------------------------------------
# Zoomed plots -- load ONLY chr1 + chr2
#
# For a large genome-wide, fine-window fst.csv, this is the important
# part: passing `chroms` here means awk filters the file down to just
# these two chromosomes BEFORE it ever reaches R, so the other 4
# chromosomes' worth of rows are never loaded into memory. Much faster
# and lighter than loading everything and filtering afterward.
# ---------------------------------------------------------
fst_data_zoom   <- load_fst(fst_csv_path, chroms = c(tbb_iso1_chrom, tbb_iso2_chrom), pair = pair)
fst_values_zoom <- fst_data_zoom$values

# sanity checks -- confirm the pair name is right and see what's available
n_distinct(fst_values_zoom$sample_pair)
sort(unique(fst_values_zoom$sample_pair))

# ---- Zoomed isotype-1 (bp scale, annotated) ----
plot_fst_zoom(
  fst_values_zoom,
  pair = pair,
  chrom = tbb_iso1_chrom,
  filter_region = tbb_iso1_region,
  gene_region = tbb_iso1_region,
  candidate_sites = tbb_iso1_sites,
  title = "Isotype-1 \u03b2-tubulin -- fine-scale FST",
  save_path = "fst_iso1_zoom.png"
)

# ---- Zoomed isotype-2 (bp scale, annotated) ----
plot_fst_zoom(
  fst_values_zoom,
  pair = pair,
  chrom = tbb_iso2_chrom,
  filter_region = tbb_iso2_region,
  gene_region = tbb_iso2_region,
  candidate_sites = tbb_iso2_sites,
  title = "Isotype-2 \u03b2-tubulin -- fine-scale FST",
  save_path = "fst_iso2_zoom.png"
)

# ---- Whole chromosome I, isotype-1 gene marked but not filtered ----
plot_fst_zoom(
  fst_values_zoom,
  pair = pair,
  chrom = tbb_iso1_chrom,
  filter_region = NULL,           # no filtering, show everything on chr1
  gene_region = tbb_iso1_region,  # still shade where the gene is
  candidate_sites = tbb_iso1_sites,
  loess_span = 0.1,               # smaller span since there's more data spread across a whole chromosome
  title = "Chromosome I -- fine-scale FST",
  save_path = "fst_chr1_whole.png"
)

# ---- Zoomed strictly to isotype-1 gene boundaries, with small flank ----
plot_fst_zoom(
  fst_values_zoom,
  pair = pair,
  chrom = tbb_iso1_chrom,
  filter_region = tbb_iso1_region,   # gene start/end
  flank = 1000,                      # +/- 1kb of context on each side
  gene_region = tbb_iso1_region,     # shaded rectangle = actual gene span
  loess_span = 0.3,
  title = "Isotype-1 \u03b2-tubulin gene boundaries -- fine-scale FST",
  save_path = "fst_iso1_gene_boundary_zoom.png"
)

fst_values_zoom %>%
  filter(sample_pair == pair,
         chrom == tbb_iso1_chrom,
         fst > 0.8) %>%
  arrange(desc(fst)) %>%
  select(chrom, start, end, fst) %>%
  print(n = 20)
