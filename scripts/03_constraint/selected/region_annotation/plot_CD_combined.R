############################################################
## GERP pie plots – C & D combined (single PDF, one legend)
## Input : gerp_constrained_master.csv
## Output: GERP_pie_C_D_combined.pdf
############################################################

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(tidyr)

## =========================================================
## 1. Region colors & order (black → white, publication-safe)
## =========================================================
region_colors <- c(
  "0-fold"      = "#000000",
  "4-fold"      = "#2F2F2F",
  "other-fold"  = "#4F4F4F",
  "Intron"      = "#7A7A7A",
  "UTR"         = "#A6A6A6",
  "UpDown5K"    = "#C8C8C8",
  "TE"          = "#E0E0E0",
  "Inter-genic" = "#F5F5F5"
)

region_levels <- names(region_colors)
thresholds    <- c("GERP>0", "GERP>2", "GERP>4", "GERP>6")

## =========================================================
## 2. Read & prepare data
## =========================================================
df <- read_csv("gerp_constrained_master.csv")

df <- df %>%
  filter(region %in% region_levels) %>%
  mutate(
    region    = factor(region, levels = region_levels),
    threshold = factor(threshold, levels = thresholds)
  )

## =========================================================
## 3. Percentage calculation
## =========================================================
df_base <- df %>%
  group_by(subgenome, threshold) %>%
  mutate(
    total   = sum(constrained),
    frac    = constrained / total,
    percent = frac * 100
  ) %>%
  ungroup()

## =========================================================
## 4. Single pie plot function
## =========================================================
plot_one_pie <- function(subg, th) {

  d <- df_base %>%
    filter(subgenome == subg, threshold == th) %>%
    complete(
      region = factor(region_levels, levels = region_levels),
      fill   = list(constrained = 0, frac = 0, percent = 0)
    ) %>%
    arrange(region) %>%
    mutate(
      ymax = cumsum(frac),
      ymin = lag(ymax, default = 0),
      mid  = (ymin + ymax) / 2
    )

  ggplot(d, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = region)) +
    geom_rect(color = "black", linewidth = 0.4) +
    coord_polar(theta = "y") +
    geom_segment(
      aes(x = 1, xend = 1.08, y = mid, yend = mid),
      color = "black",
      linewidth = 0.35
    ) +
    geom_text(
      aes(x = 1.15, y = mid, label = paste0(round(percent), "%")),
      size = 3
    ) +
    scale_fill_manual(values = region_colors, drop = FALSE) +
    labs(title = th) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 11, hjust = 0.5),
      legend.position = "none"
    )
}

## =========================================================
## 5. Build C / D columns
## =========================================================
plot_column <- function(subg, col_title) {

  pies <- lapply(thresholds, function(th) plot_one_pie(subg, th))

  wrap_plots(pies, ncol = 1) +
    plot_annotation(title = col_title) &
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
    )
}

p_C <- plot_column("C", "C subgenome")
p_D <- plot_column("D", "D subgenome")

## =========================================================
## 6. Combine C & D, collect legend (top-right)
## =========================================================
final_plot <- (p_C | p_D) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title    = element_blank()
  )

## =========================================================
## 7. Save final figure
## =========================================================
ggsave(
  "GERP_pie_C_D_combined.pdf",
  final_plot,
  width = 11,
  height = 10
)

final_plot
