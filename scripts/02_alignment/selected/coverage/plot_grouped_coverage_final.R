#!/usr/bin/env Rscript
# Observed coverage by subfamily for one reference-coordinate alignment.
# Usage: Rscript plot_grouped_coverage_final.R CLASS_LIST COVERAGE_TSV OUTPUT_PDF [TITLE]
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: Rscript plot_grouped_coverage_final.R CLASS_LIST COVERAGE_TSV OUTPUT_PDF [TITLE]")
}
title <- if (length(args) == 4) args[[4]] else "Observed coverage by subfamily"

class_data <- read.table(args[[1]], header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("Species", "Subfamily"))
cov_data <- read.delim(args[[2]], header = TRUE, stringsAsFactors = FALSE)
if (!all(c("Species", "ObsCoverage") %in% names(cov_data))) {
  stop("Coverage table needs Species and ObsCoverage columns")
}
if (anyDuplicated(class_data$Species) || anyDuplicated(cov_data$Species)) {
  stop("Species identifiers must be unique in both input tables")
}
if (!setequal(class_data$Species, cov_data$Species)) {
  stop("Classification and coverage species sets differ")
}

plot_data <- class_data %>%
  left_join(cov_data, by = "Species") %>%
  mutate(Species = factor(Species, levels = class_data$Species),
         ObsCoverage = as.numeric(ObsCoverage))
if (anyNA(plot_data$ObsCoverage) || any(plot_data$ObsCoverage < 0 | plot_data$ObsCoverage > 100)) {
  stop("Observed coverage must be numeric percentages between 0 and 100")
}

subfamily_colors <- c(
  "Bambusoideae" = "#1f77b4", "Oryzoideae" = "#ff7f0e",
  "Pooideae" = "#2ca02c", "Arundinoideae" = "#d62728",
  "Chloridoideae" = "#9467bd", "Panicoideae" = "#8c564b",
  "Anomochlooideae" = "#e377c2", "Pharoideae" = "#7f7f7f"
)
if (any(!plot_data$Subfamily %in% names(subfamily_colors))) {
  stop("Unknown subfamily in classification table")
}

# Keep the supplied 0-80% axis when possible, and expand it if needed.
ymax <- max(80, ceiling(max(plot_data$ObsCoverage) / 20) * 20)
p <- ggplot(plot_data, aes(x = Species, y = ObsCoverage, fill = Subfamily)) +
  geom_col(width = 0.7, color = NA) +
  scale_fill_manual(values = subfamily_colors) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax, by = 20),
                     labels = number_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(title = title, x = NULL, y = "Observed coverage (%)") +
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               size = 5, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey", linewidth = 0.3),
    legend.margin = margin(4, 4, 4, 4),
    legend.key.size = grid::unit(3, "mm"),
    legend.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(ncol = 2))

ggsave(args[[3]], plot = p, width = 180, height = 100, units = "mm", dpi = 600)
