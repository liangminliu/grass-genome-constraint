#!/usr/bin/env Rscript
# GERP >2 region composition for the three focal grass species.
# Input: headerless TSV with species, region, count columns.
# Default input: all_gt2_CDS.tsv
# Default output: GERP_gt2_species_pies.pdf
# Usage: Rscript plot_gt2_multi_species.R [input.tsv] [output.pdf]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2) stop("Usage: Rscript plot_gt2_multi_species.R [input.tsv] [output.pdf]")
input <- if (length(args) >= 1) args[[1]] else "all_gt2_CDS.tsv"
output <- if (length(args) >= 2) args[[2]] else "GERP_gt2_species_pies.pdf"

species_order <- c("P.edulis", "Z.may", "A.tauschii")
region_levels <- c("CDS", "Intron", "UTR", "UpDown5K", "Inter-genic")
region_colors <- c(
  "CDS" = "#2F2F2F",
  "Intron" = "#7A7A7A",
  "UTR" = "#A6A6A6",
  "UpDown5K" = "#C8C8C8",
  "Inter-genic" = "#F0F0F0"
)

raw <- read_tsv(
  input,
  col_names = c("species", "region", "count"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)
if (ncol(raw) != 3L || nrow(raw) == 0L) stop("Input must contain three columns and at least one row")
if (anyNA(raw$species) || anyNA(raw$region) || anyNA(raw$count)) stop("Input contains missing values")
if (any(!raw$species %in% species_order)) stop("Unknown species identifier in input")
if (any(!raw$region %in% region_levels)) stop("Unknown region category in input")
if (anyDuplicated(raw[c("species", "region")])) stop("Duplicate species-region pair in input")

counts <- as.numeric(gsub(",", "", raw$count, fixed = TRUE))
if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) {
  stop("Counts must be non-negative integers; thousands separators are allowed")
}

df <- raw %>%
  mutate(
    count = counts,
    species = factor(species, levels = species_order),
    region = factor(region, levels = region_levels)
  ) %>%
  complete(species, region, fill = list(count = 0)) %>%
  group_by(species) %>%
  mutate(total = sum(count), frac = if_else(total > 0, count / total, 0),
         percent = 100 * frac) %>%
  ungroup()

if (any(df$total == 0)) stop("Each focal species needs a positive total count")

plot_one_species <- function(sp) {
  d <- df %>%
    filter(species == sp) %>%
    arrange(region) %>%
    mutate(
      ymax = cumsum(frac),
      ymin = lag(ymax, default = 0),
      mid = (ymin + ymax) / 2,
      label = if_else(count > 0, paste0(round(percent), "%"), "")
    )

  ggplot(d, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = region)) +
    geom_rect(color = "black", linewidth = 0.4) +
    geom_segment(aes(x = 1, xend = 1.06, y = mid, yend = mid),
                 linewidth = 0.35) +
    geom_text(aes(x = 1.12, y = mid, label = label), size = 3) +
    coord_polar(theta = "y", clip = "off") +
    scale_fill_manual(values = region_colors, breaks = region_levels, drop = FALSE) +
    labs(title = sp, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"))
}

plots <- lapply(species_order, plot_one_species)
final_plot <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.title = element_blank())

ggsave(output, final_plot, width = 12, height = 4)
print(final_plot)
