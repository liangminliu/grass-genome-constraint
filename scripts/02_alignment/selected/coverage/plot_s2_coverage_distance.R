#!/usr/bin/env Rscript

# ================================================================
# Coverage–distance plots | Final V3 (independent panel scales)
# References: P. edulis, Z. mays, A. tauschii
# Outputs: Coverage d–f, ObsCoverage g–i
# Analysis: all valid species; axes never change the fitted model.
# ================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
  library(grid)
})

# -------------------- 1. Global parameters ---------------------
OUTDIR <- "coverage_tripanel_final_v3"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUTGROUPS <- c("Acom", "Jeff", "Ssto")
FIG_WIDTH  <- 14.4  # inches, for 3 panels
FIG_HEIGHT <- 4.65
POINT_SIZE <- 1.75
POINT_ALPHA <- 0.82
POINT_STROKE <- 0.22
LINE_WIDTH <- 0.75
FONT_FAMILY <- "sans"       # portable; use Arial if available

# Manual ranges are INITIAL suggestions, not a species filter.
# 'manual_expand' retains ALL observations by widening boundaries
# whenever a point falls outside a user-specified range.
AXIS_MODE <- "manual_expand"

# Optional black transparent silhouettes: PNG or SVG.
# Keep empty unless you have actual local image files.
SILHOUETTE_FILES <- c(ped = "", Zmay = "", Atau = "")
SHOW_SILHOUETTES <- FALSE

# ---------------- 2. Reference configurations ------------------
cfgs <- list(
  list(
    ref_id = "ped", ref_label = "Phyllostachys edulis",
    color = "#3B6FB6",
    roast_dir = "97sp_ped.results/roast",
    tree_file = "97sp_ped.results/tree/neutral.tre",
    letter_cov = "A", letter_obs = "A"
  ),
  list(
    ref_id = "Zmay", ref_label = "Zea mays",
    color = "#C17A3A",
    roast_dir = "../Zmay/Zmay_v5_results/roast",
    tree_file = "../Zmay/Zmay_v5_results/tree/neutral.tre",
    letter_cov = "B", letter_obs = "B"
  ),
  list(
    ref_id = "Atau", ref_label = "Aegilops tauschii",
    color = "#3F8664",
    roast_dir = "../Atau_msa/results/roast",
    tree_file = "../Atau_msa/results/tree/neutral.tre",
    letter_cov = "C", letter_obs = "C"
  )
)

# ---------------- 3. INDIVIDUAL AXIS SETTINGS -------------------
# Change the x/y ranges here, separately for each species and metric.
# Range limits are extended automatically if actual points exceed them.
# tick_x/tick_y are major tick spacing, in original data units.
AXIS_SETTINGS <- list(
  ped = list(
    Coverage    = list(x = c(0.16, 1.00), y = c(0.00, 0.24),
                       tick_x = 0.2, tick_y = 0.05),
    ObsCoverage = list(x = c(0.16, 1.00), y = c(0.10, 0.80),
                       tick_x = 0.2, tick_y = 0.10)
  ),
  Zmay = list(
    Coverage    = list(x = c(0.20, 1.10), y = c(0.00, 0.19),
                       tick_x = 0.2, tick_y = 0.05),
    ObsCoverage = list(x = c(0.20, 1.10), y = c(0.05, 0.66),
                       tick_x = 0.2, tick_y = 0.10)
  ),
  Atau = list(
    Coverage    = list(x = c(0.05, 1.25), y = c(0.00, 0.68),
                       tick_x = 0.2, tick_y = 0.10),
    ObsCoverage = list(x = c(0.05, 1.25), y = c(0.00, 0.45),
                       tick_x = 0.2, tick_y = 0.10)
  )
)

# --------------------- 4. Input parsing -------------------------
parse_cov_file <- function(f, reference_id) {
  x <- readLines(f, warn = FALSE)
  # Overall reference row: ped*, Zmay*, Atau*, etc.
  pattern <- paste0("^[[:space:]]*", reference_id, "\\*[[:space:]]+")
  hits <- grep(pattern, x)
  if (!length(hits)) return(NULL)

  fields <- strsplit(trimws(x[hits[1]]), "[[:space:]]+")[[1]]
  if (length(fields) < 6) {
    warning("Fewer than six columns in ", f)
    return(NULL)
  }
  values <- suppressWarnings(as.numeric(fields[2:6]))
  if (any(!is.finite(values))) {
    warning("Non-numeric coverage data in ", f)
    return(NULL)
  }

  tibble(
    Species = sub("\\.tsv$", "", basename(f), ignore.case = TRUE),
    SourceLength = values[1],
    ObsLength = values[2],
    AlignedPos = values[3],
    Coverage = values[4],
    ObsCoverage = values[5]
  )
}

# ------------------ 5. Pearson + OLS on ALL data ---------------
get_stats <- function(df, metric, ref_id) {
  x <- df$Distance
  y <- df[[metric]]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || sd(x) <= 0 || sd(y) <= 0) {
    stop("Insufficient variation: ", ref_id, ", ", metric)
  }
  fit <- lm(y ~ x)
  pear <- cor.test(x, y, method = "pearson")
  spear <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(
    Reference = ref_id, metric = metric, n = length(x),
    pearson_r = unname(pear$estimate), pearson_p = pear$p.value,
    spearman_r = unname(spear$estimate), spearman_p = spear$p.value,
    slope = unname(coef(fit)[2]),
    intercept = unname(coef(fit)[1]),
    r_squared = summary(fit)$r.squared
  )
}

# -------------------- 6. Load one reference ---------------------
load_reference <- function(cfg) {
  class_file <- file.path(cfg$roast_dir, "94sp.pop.list")
  cov_dir <- file.path(cfg$roast_dir, "mafPairCoverage_results")
  for (path in c(class_file, cfg$tree_file)) {
    if (!file.exists(path)) stop("Missing input file: ", path)
  }
  if (!dir.exists(cov_dir)) stop("Missing coverage directory: ", cov_dir)

  class_data <- read.table(
    class_file, header = FALSE, stringsAsFactors = FALSE,
    comment.char = "#", fill = TRUE
  )
  if (ncol(class_data) < 2) stop("Expected 2 columns in ", class_file)
  class_data <- tibble(
    Species = as.character(class_data[[1]]),
    Subfamily = as.character(class_data[[2]])
  )
  if (anyDuplicated(class_data$Species)) {
    stop("Duplicate species in ", class_file)
  }

  fs <- list.files(cov_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(fs)) stop("No .tsv coverage files in ", cov_dir)
  cov_data <- bind_rows(lapply(fs, parse_cov_file, reference_id = cfg$ref_id))
  if (!nrow(cov_data)) stop("No coverage parsed: ", cfg$ref_id)
  if (anyDuplicated(cov_data$Species)) {
    stop("Duplicate species from coverage file names: ", cfg$ref_id)
  }

  tr <- ape::read.tree(cfg$tree_file)
  if (!(cfg$ref_id %in% tr$tip.label)) {
    stop("Reference absent from tree: ", cfg$ref_id)
  }
  dist_matrix <- ape::cophenetic.phylo(tr)
  dist_data <- tibble(
    Species = colnames(dist_matrix),
    Distance = as.numeric(dist_matrix[cfg$ref_id, ])
  )

  dat <- cov_data %>%
    inner_join(class_data, by = "Species") %>%
    inner_join(dist_data, by = "Species") %>%
    filter(
      Species != cfg$ref_id,
      !(Species %in% OUTGROUPS),
      tolower(Subfamily) != "outgroup",
      is.finite(Distance),
      is.finite(Coverage), is.finite(ObsCoverage),
      between(Coverage, 0, 1),
      between(ObsCoverage, 0, 1)
    ) %>%
    arrange(Distance) %>%
    mutate(Reference = cfg$ref_id, .before = 1)

  if (nrow(dat) < 3) stop("Too few valid species: ", cfg$ref_id)
  if (anyDuplicated(dat$Species)) stop("Duplicate merged species: ", cfg$ref_id)

  stats <- bind_rows(
    get_stats(dat, "Coverage", cfg$ref_id),
    get_stats(dat, "ObsCoverage", cfg$ref_id)
  )
  message("[DATA] ", cfg$ref_id, ": ", nrow(dat), " species; ",
          "Coverage max=", round(max(dat$Coverage), 3),
          "; ObsCoverage max=", round(max(dat$ObsCoverage), 3))
  list(data = dat, stats = stats, cfg = cfg)
}

# -------------------- 7. Safe axis limits -----------------------
resolve_axes <- function(dat, ref_id, metric) {
  ax <- AXIS_SETTINGS[[ref_id]][[metric]]
  if (is.null(ax)) stop("Missing axis definition: ", ref_id, "/", metric)
  xlim <- as.numeric(ax$x)
  ylim <- as.numeric(ax$y)
  if (length(xlim) != 2 || length(ylim) != 2 ||
      any(!is.finite(c(xlim, ylim))) ||
      diff(xlim) <= 0 || diff(ylim) <= 0 ||
      ax$tick_x <= 0 || ax$tick_y <= 0) {
    stop("Invalid axis definition: ", ref_id, "/", metric)
  }
  if (AXIS_MODE != "manual_expand") stop("Unsupported AXIS_MODE")

  # Preserve ALL original datapoints in the graph, including outliers.
  xx <- range(dat$Distance)
  yy <- range(dat[[metric]])
  expanded <- xx[1] < xlim[1] || xx[2] > xlim[2] ||
              yy[1] < ylim[1] || yy[2] > ylim[2]
  if (expanded) {
    pad_x <- diff(xx) * 0.025
    pad_y <- diff(yy) * 0.025
    xlim <- c(min(xlim[1], xx[1] - pad_x),
              max(xlim[2], xx[2] + pad_x))
    ylim <- c(min(ylim[1], yy[1] - pad_y),
              max(ylim[2], yy[2] + pad_y))
    # Fractions and branch distances are nonnegative.
    xlim[1] <- max(0, xlim[1])
    ylim[1] <- max(0, ylim[1])
    message("[AXIS EXPANDED] ", ref_id, " / ", metric,
            "  X=[", paste(round(xlim, 3), collapse = ", "),
            "] Y=[", paste(round(ylim, 3), collapse = ", "), "]")
  }
  list(
    x = xlim, y = ylim,
    tick_x = ax$tick_x, tick_y = ax$tick_y,
    expanded = expanded
  )
}

make_breaks <- function(limits, step) {
  a <- ceiling((limits[1] - 1e-9) / step) * step
  b <- floor((limits[2] + 1e-9) / step) * step
  if (a > b) return(numeric(0))
  round(seq(a, b, by = step), 6)
}

format_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 1e-300) return("P < 1e-300")
  sprintf("P = %.2e", p)
}

# ----------- 8. Optional SVG/PNG silhouette reader -------------
read_silhouette <- function(ref_id) {
  if (!SHOW_SILHOUETTES) return(NULL)
  p <- unname(SILHOUETTE_FILES[ref_id])
  if (length(p) != 1 || is.na(p) || !nzchar(p)) return(NULL)
  if (!file.exists(p)) {
    warning("Silhouette unavailable: ", p)
    return(NULL)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    warning("Install R package 'png' to show silhouettes")
    return(NULL)
  }
  ext <- tolower(tools::file_ext(p))
  if (ext == "png") return(png::readPNG(p))
  if (ext == "svg") {
    if (!requireNamespace("rsvg", quietly = TRUE)) {
      warning("Install R package 'rsvg' to use SVG silhouettes")
      return(NULL)
    }
    return(png::readPNG(rsvg::rsvg_png(p, width = 400, height = 400)))
  }
  warning("Silhouette must be PNG or SVG: ", p)
  NULL
}

# -------------------- 9. Scientific figure ----------------------
make_panel <- function(rr, metric, panel_letter, panel_i) {
  dat <- rr$data
  cfg <- rr$cfg
  st <- rr$stats %>% filter(.data$metric == !!metric)
  ax <- resolve_axes(dat, cfg$ref_id, metric)
  xl <- ax$x
  yl <- ax$y
  dx <- diff(xl)
  dy <- diff(yl)

  # Upper-right annotation in relative coordinates, so each axis
  # can be adjusted independently without disrupting typography.
  annot_x <- xl[2] - 0.038 * dx
  name_y  <- yl[2] - 0.045 * dy
  r_y     <- yl[2] - 0.135 * dy
  p_y     <- yl[2] - 0.215 * dy

  y_axis_label <- if (metric == "Coverage") {
    "Alignment coverage relative\nto reference genome"
  } else {
    "Observed-base alignment\ncoverage"
  }

  p <- ggplot(dat, aes(x = Distance, y = .data[[metric]])) +
    # OLS coefficients were computed with ALL valid observations.
    # coord_cartesian clip='on' prevents the line leaving the frame.
    geom_abline(
      intercept = st$intercept[[1]], slope = st$slope[[1]],
      color = cfg$color, linewidth = LINE_WIDTH,
      alpha = 0.98, show.legend = FALSE
    ) +
    geom_point(
      shape = 21, fill = cfg$color, color = "#202020",
      size = POINT_SIZE, stroke = POINT_STROKE,
      alpha = POINT_ALPHA, show.legend = FALSE
    ) +
    annotate(
      "text", x = annot_x, y = name_y,
      label = cfg$ref_label, family = FONT_FAMILY,
      hjust = 1, vjust = 1, size = 4.25,
      fontface = "bold.italic", color = "black"
    ) +
    annotate(
      "text", x = annot_x, y = r_y,
      label = sprintf("Pearson r = %.2f", st$pearson_r[[1]]),
      family = FONT_FAMILY,
      hjust = 1, vjust = 1, size = 3.55, color = "black"
    ) +
    annotate(
      "text", x = annot_x, y = p_y,
      label = format_p(st$pearson_p[[1]]),
      family = FONT_FAMILY,
      hjust = 1, vjust = 1, size = 3.55, color = "black"
    ) +
    scale_x_continuous(
      breaks = make_breaks(xl, ax$tick_x),
      labels = function(x) sprintf("%.1f", x),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = make_breaks(yl, ax$tick_y),
      labels = function(y) sprintf("%.2f", y),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = xl, ylim = yl, clip = "on") +
    labs(
      tag = panel_letter,
      x = if (panel_i == 2) "Neutral branch length distance to reference" else NULL,
      y = if (panel_i == 1) y_axis_label else NULL
    ) +
    theme_classic(base_size = 11, base_family = FONT_FAMILY) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.65),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = grid::unit(2.0, "mm"),
      axis.text = element_text(color = "black", size = 9.5),
      axis.title = element_text(color = "black", size = 11.0),
      axis.title.x = element_text(margin = margin(t = 9)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.tag = element_text(face = "bold", size = 19,
                              family = FONT_FAMILY, color = "black"),
      plot.tag.position = c(0.005, 1.0),
      plot.margin = margin(t = 13, r = 5, b = 5, l = 7)
    )

  icon <- read_silhouette(cfg$ref_id)
  if (!is.null(icon)) {
    # Center-left placement; avoids the labels at upper right.
    # For data overlapping this area, disable silhouettes.
    p <- p + annotation_custom(
      grob = grid::rasterGrob(icon, interpolate = TRUE),
      xmin = xl[1] + 0.28 * dx, xmax = xl[1] + 0.37 * dx,
      ymin = yl[2] - 0.19 * dy, ymax = yl[2] - 0.02 * dy
    )
  }

  # Extra sanity check: NO observed point is outside plotted range.
  outside <- sum(dat$Distance < xl[1] | dat$Distance > xl[2] |
                 dat[[metric]] < yl[1] | dat[[metric]] > yl[2])
  if (outside != 0L) stop("Some points were clipped in ", cfg$ref_id, "/", metric)

  audit <- tibble(
    Reference = cfg$ref_id, metric = metric,
    panel = panel_letter, n = nrow(dat),
    x_min = xl[1], x_max = xl[2],
    y_min = yl[1], y_max = yl[2],
    adjusted_to_show_all = ax$expanded,
    n_points_outside_axes = outside
  )
  list(plot = p, audit = audit)
}

# ---------------------- 10. Figure output -----------------------
results <- lapply(cfgs, load_reference)
all_data <- bind_rows(lapply(results, `[[`, "data"))
all_stats <- bind_rows(lapply(results, `[[`, "stats"))
readr::write_tsv(all_data,
  file.path(OUTDIR, "all_species_coverage_distance.tsv"))
readr::write_tsv(all_stats,
  file.path(OUTDIR, "all_species_statistics.tsv"))
print(all_stats, n = Inf)

make_tripanel <- function(metric) {
  assembled <- lapply(seq_along(results), function(i) {
    rr <- results[[i]]
    letter <- if (metric == "Coverage") rr$cfg$letter_cov else rr$cfg$letter_obs
    make_panel(rr, metric, letter, i)
  })
  fig <- patchwork::wrap_plots(
    lapply(assembled, `[[`, "plot"), ncol = 3
  ) + patchwork::plot_layout(widths = c(1, 1, 1))
  audits <- bind_rows(lapply(assembled, `[[`, "audit"))
  list(fig = fig, audits = audits)
}

save_tripanel <- function(metric) {
  obj <- make_tripanel(metric)
  basename_out <- paste0(metric, "_tripanel")
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    filename = file.path(OUTDIR, paste0(basename_out, ".pdf")),
    plot = obj$fig, width = FIG_WIDTH, height = FIG_HEIGHT,
    units = "in", device = pdf_device, bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(OUTDIR, paste0(basename_out, ".png")),
    plot = obj$fig, width = FIG_WIDTH, height = FIG_HEIGHT,
    units = "in", dpi = 600, bg = "white"
  )
  message("[SAVED] ", basename_out, ".pdf / .png")
  obj$audits
}

audit_cov <- save_tripanel("Coverage")
audit_obs <- save_tripanel("ObsCoverage")
axis_audit <- bind_rows(audit_cov, audit_obs)
readr::write_tsv(axis_audit, file.path(OUTDIR, "axis_ranges_used.tsv"))
print(axis_audit, n = Inf)
message("[DONE] Output directory: ", normalizePath(OUTDIR))
