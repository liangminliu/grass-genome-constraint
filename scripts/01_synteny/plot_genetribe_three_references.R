#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

summary_file <- ifelse(length(args) >= 1, args[1], "genetribe_main_summary.tsv")
order_file   <- ifelse(length(args) >= 2, args[2], "draw.sample.list")
outdir       <- ifelse(length(args) >= 3, args[3], "genetribe_three_reference_plots")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Main switches
# =========================================================
# TRUE: self-reference values are not drawn.
#       But the 3 ref slots are still kept, so each species keeps 3 bar positions.
REMOVE_SELF_VALUES <- TRUE

# TRUE: keep every species in draw.sample.list / tree order,
#       including species without OK results, such as Coix_aquatica.
KEEP_TREE_EMPTY_SPECIES <- TRUE

# Horizontal bar geometry with manual y coordinates
REF_BAR_STEP <- 0.22
BAR_HEIGHT   <- 0.22

SYNTENIC_XMAX    <- 0.95
COMPOSITION_XMAX <- 1.0

# =========================================================
# Helper functions
# =========================================================
norm_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  tolower(x)
}

std_species_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x
}

abbr_species <- function(x) {
  x <- std_species_name(x)
  parts <- unlist(strsplit(x, "_"))
  parts <- parts[parts != ""]
  if (length(parts) >= 2) {
    paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
  } else {
    gsub("_", " ", x)
  }
}

blank_labeler <- function(x) rep("", length(x))

pct_int <- function(x) paste0(round(x * 100), "%")

make_breaks_syntenic <- function() {
  c(0, 0.2, 0.4, 0.6, 0.8, 0.95)
}

make_breaks_composition <- function() {
  c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
}

first_non_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  y <- x[!is.na(x)]
  if (length(y) == 0) NA_real_ else y[1]
}

# =========================================================
# Output size
# =========================================================
WIDTH_LABEL_SINGLE   <- 3.55
WIDTH_NOLABEL_SINGLE <- 2.45
HEIGHT_SINGLE        <- 11.40

WIDTH_LABEL_COMBO    <- 6.30
WIDTH_NOLABEL_COMBO  <- 4.70
HEIGHT_COMBO         <- 11.40

# =========================================================
# Read data
# =========================================================
dat <- read.delim(
  summary_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tree_order_raw <- readLines(order_file, warn = FALSE)
tree_order_raw <- std_species_name(tree_order_raw)
tree_order_raw <- tree_order_raw[nzchar(tree_order_raw)]
tree_order_raw <- tree_order_raw[!grepl("^#", tree_order_raw)]

if (!"result_source_mode" %in% colnames(dat)) {
  dat$result_source_mode <- NA_character_
}

if (!"status" %in% colnames(dat)) {
  stop("Input summary table must contain a 'status' column.")
}

dat <- dat %>%
  mutate(
    query_species = std_species_name(query_species),
    ref_species   = std_species_name(ref_species),
    query_key     = norm_key(query_species),
    ref_key       = norm_key(ref_species)
  )

dat_ok <- dat %>%
  filter(grepl("^OK", status))

num_cols <- c(
  "prop_syntenic_genes_ref",
  "prop_syntenic_genes_query",
  "prop_singleton_ref",
  "prop_singleton_query"
)

for (cc in intersect(num_cols, colnames(dat_ok))) {
  dat_ok[[cc]] <- suppressWarnings(as.numeric(dat_ok[[cc]]))
}

dat_ok <- dat_ok %>%
  mutate(
    is_self = query_key == ref_key |
      result_source_mode == "synthetic_self" |
      grepl("^OK_self", status)
  )

if (REMOVE_SELF_VALUES) {
  dat_ok <- dat_ok %>% filter(!is_self)
}

# =========================================================
# Species order from tree / draw.sample.list
# Important:
# 1. Keep all species in tree_order_raw.
# 2. Therefore Coix_aquatica will keep its tree position even if no OK result exists.
# 3. Case mismatch and space/underscore mismatch are handled by norm_key().
# =========================================================
order_df <- data.frame(
  order_in_tree = seq_along(tree_order_raw),
  query_species_tree = tree_order_raw,
  query_key = norm_key(tree_order_raw),
  stringsAsFactors = FALSE
) %>%
  distinct(query_key, .keep_all = TRUE)

query_lookup <- dat %>%
  distinct(query_key, query_species) %>%
  group_by(query_key) %>%
  summarise(
    query_species_data = first(query_species),
    .groups = "drop"
  )

species_order_df <- order_df %>%
  left_join(query_lookup, by = "query_key") %>%
  mutate(
    query_species_plot = ifelse(
      !is.na(query_species_data),
      query_species_tree,
      query_species_tree
    )
  )

# Add extra query species that are present in OK results but absent from tree list.
# They are appended to the bottom.
extra_queries_df <- dat_ok %>%
  distinct(query_key, query_species) %>%
  anti_join(species_order_df %>% distinct(query_key), by = "query_key") %>%
  arrange(query_species)

if (nrow(extra_queries_df) > 0) {
  extra_df <- data.frame(
    order_in_tree = max(species_order_df$order_in_tree, na.rm = TRUE) + seq_len(nrow(extra_queries_df)),
    query_species_tree = extra_queries_df$query_species,
    query_key = extra_queries_df$query_key,
    query_species_data = extra_queries_df$query_species,
    query_species_plot = extra_queries_df$query_species,
    stringsAsFactors = FALSE
  )
  species_order_df <- bind_rows(species_order_df, extra_df)
}

species_order_df <- species_order_df %>%
  arrange(order_in_tree) %>%
  distinct(query_key, .keep_all = TRUE)

# Manual y center:
# First species in tree is plotted at the top.
n_species <- nrow(species_order_df)

species_meta <- species_order_df %>%
  mutate(
    y_center = n_species - row_number() + 1,
    y_label  = sapply(query_species_plot, abbr_species, USE.NAMES = FALSE)
  )

# =========================================================
# Reference order
# Keep exactly 3 reference slots for each species.
# =========================================================
target_refs_raw <- c(
  "Phyllostachys_edulis",
  "Zea_mays",
  "Aegilops_tauschii"
)

ref_title_map <- c(
  "Phyllostachys_edulis" = "P. edulis",
  "Zea_mays"             = "Z. mays",
  "Aegilops_tauschii"    = "A. tauschii"
)

ref_display <- setNames(
  sapply(
    target_refs_raw,
    function(rr) {
      if (rr %in% names(ref_title_map)) {
        ref_title_map[[rr]]
      } else {
        abbr_species(rr)
      }
    },
    USE.NAMES = FALSE
  ),
  target_refs_raw
)

ref_meta <- data.frame(
  ref_species = target_refs_raw,
  ref_key     = norm_key(target_refs_raw),
  ref_index   = seq_along(target_refs_raw),
  stringsAsFactors = FALSE
) %>%
  mutate(
    ref_label = unname(ref_display[ref_species]),
    ref_offset = ((length(target_refs_raw) + 1) / 2 - ref_index) * REF_BAR_STEP
  )

# Full plotting grid:
# Every species in tree × every reference.
# Missing values remain NA and are not drawn, but the vertical space is retained.
row_grid <- tidyr::expand_grid(
  species_meta %>%
    select(
      order_in_tree,
      query_species_plot,
      query_key,
      y_center,
      y_label
    ),
  ref_meta
) %>%
  mutate(
    y = y_center + ref_offset,
    ref_label = factor(ref_label, levels = unname(ref_display[target_refs_raw]))
  )

# =========================================================
# Palettes
# =========================================================
ref_palette <- c(
  "P. edulis"   = "#B7CBE2",
  "Z. mays"     = "#D8C1B1",
  "A. tauschii" = "#B8CDAF"
)

comp_palette <- c(
  "Syntenic"  = "#8FAFD0",
  "Singleton" = "#CDBA8F",
  "Other"     = "#E6E6E6"
)

# =========================================================
# Theme
# =========================================================
base_theme <- theme_minimal(base_size = 8.6) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      colour = "#F0F0F0",
      linewidth = 0.25,
      linetype = "solid"
    ),

    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8.6, colour = "#333333", margin = margin(t = 4)),

    axis.text.x = element_text(size = 7.4, colour = "#333333"),
    axis.text.y = element_text(size = 5.9, colour = "#333333", face = "italic"),

    axis.ticks.x = element_line(colour = "#888888", linewidth = 0.25),
    axis.ticks.y = element_line(colour = "#888888", linewidth = 0.25),
    axis.ticks.length = grid::unit(1.8, "pt"),

    axis.line.x = element_line(colour = "#888888", linewidth = 0.25),
    axis.line.y = element_line(colour = "#888888", linewidth = 0.25),

    panel.border = element_rect(fill = NA, colour = "#D0D0D0", linewidth = 0.28),

    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = element_text(size = 8.0, colour = "#333333"),
    legend.text = element_text(size = 7.5, colour = "#333333"),
    legend.key.height = grid::unit(3.5, "pt"),
    legend.key.width  = grid::unit(10, "pt"),

    plot.title = element_text(size = 10.6, face = "bold", colour = "#222222"),
    plot.subtitle = element_text(size = 7.8, colour = "#555555"),
    plot.caption = element_text(size = 6.8, colour = "#555555", hjust = 0),

    plot.margin = margin(4, 4, 4, 4)
  )

# =========================================================
# Build plotting tables
# =========================================================
make_metric_df <- function(df, metric_col) {
  raw <- df %>%
    transmute(
      query_key = query_key,
      ref_key   = ref_key,
      value     = .data[[metric_col]]
    ) %>%
    group_by(query_key, ref_key) %>%
    summarise(
      value = first_non_na(value),
      .groups = "drop"
    )

  out <- row_grid %>%
    left_join(raw, by = c("query_key", "ref_key")) %>%
    mutate(
      ref_label = factor(ref_label, levels = unname(ref_display[target_refs_raw]))
    )

  out
}

make_comp_df <- function(df, syn_col, sing_col) {
  raw <- df %>%
    transmute(
      query_key = query_key,
      ref_key   = ref_key,
      Syntenic  = .data[[syn_col]],
      Singleton = .data[[sing_col]]
    ) %>%
    group_by(query_key, ref_key) %>%
    summarise(
      Syntenic  = first_non_na(Syntenic),
      Singleton = first_non_na(Singleton),
      .groups = "drop"
    ) %>%
    mutate(
      Other = ifelse(
        is.na(Syntenic) | is.na(Singleton),
        NA_real_,
        pmax(0, 1 - Syntenic - Singleton)
      )
    )

  out <- row_grid %>%
    left_join(raw, by = c("query_key", "ref_key")) %>%
    mutate(
      ref_label = factor(ref_label, levels = unname(ref_display[target_refs_raw]))
    ) %>%
    pivot_longer(
      cols = c(Syntenic, Singleton, Other),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component = factor(component, levels = c("Syntenic", "Singleton", "Other"))
    )

  out
}

# =========================================================
# Plot functions
# =========================================================
make_metric_plot <- function(metric_df,
                             title_text,
                             x_title,
                             show_y = TRUE,
                             xlim_max = SYNTENIC_XMAX) {

  y_breaks <- species_meta$y_center
  y_labels <- if (show_y) species_meta$y_label else rep("", nrow(species_meta))
  y_limits <- c(min(species_meta$y_center) - 0.55, max(species_meta$y_center) + 0.55)

  p <- ggplot(metric_df, aes(x = value, y = y, fill = ref_label)) +
    geom_col(
      width = BAR_HEIGHT,
      orientation = "y",
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = ref_palette,
      breaks = unname(ref_display[target_refs_raw]),
      name = "Reference",
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, xlim_max),
      breaks = make_breaks_syntenic(),
      labels = pct_int,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = title_text,
      x = x_title,
      y = NULL
    ) +
    base_theme +
    guides(
      fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)
    ) +
    theme(
      legend.position = "top",
      axis.text.y = if (show_y) {
        element_text(size = 5.9, colour = "#333333", face = "italic")
      } else {
        element_blank()
      }
    )

  p
}

make_comp_plot <- function(comp_df,
                           title_text,
                           x_title,
                           show_y = TRUE,
                           xlim_max = COMPOSITION_XMAX) {

  y_breaks <- species_meta$y_center
  y_labels <- if (show_y) species_meta$y_label else rep("", nrow(species_meta))
  y_limits <- c(min(species_meta$y_center) - 0.55, max(species_meta$y_center) + 0.55)

  marker_df <- comp_df %>%
    group_by(query_key, ref_key, y, ref_label) %>%
    summarise(
      has_value = any(!is.na(value)),
      .groups = "drop"
    ) %>%
    filter(has_value)

  p <- ggplot(comp_df, aes(x = value, y = y, fill = component)) +
    geom_col(
      width = BAR_HEIGHT,
      orientation = "y",
      na.rm = TRUE
    ) +
    geom_point(
      data = marker_df,
      aes(x = -0.035, y = y, colour = ref_label),
      inherit.aes = FALSE,
      shape = 15,
      size = 1.6
    ) +
    scale_fill_manual(
      values = comp_palette,
      name = "Component",
      drop = FALSE
    ) +
    scale_colour_manual(
      values = ref_palette,
      breaks = unname(ref_display[target_refs_raw]),
      name = "Reference",
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(-0.05, xlim_max),
      breaks = make_breaks_composition(),
      labels = function(x) ifelse(x < 0, "", pct_int(x)),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_text,
      x = x_title,
      y = NULL
    ) +
    base_theme +
    guides(
      colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
      fill   = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    ) +
    theme(
      legend.position = "top",
      axis.text.y = if (show_y) {
        element_text(size = 5.9, colour = "#333333", face = "italic")
      } else {
        element_blank()
      }
    )

  p
}

# =========================================================
# Prepare plot data
# =========================================================
metric_ref_df   <- make_metric_df(dat_ok, "prop_syntenic_genes_ref")
metric_query_df <- make_metric_df(dat_ok, "prop_syntenic_genes_query")

comp_ref_df <- make_comp_df(
  dat_ok,
  syn_col  = "prop_syntenic_genes_ref",
  sing_col = "prop_singleton_ref"
)

comp_query_df <- make_comp_df(
  dat_ok,
  syn_col  = "prop_syntenic_genes_query",
  sing_col = "prop_singleton_query"
)

# =========================================================
# Label plots
# =========================================================
p_syn_genome_label <- make_metric_plot(
  metric_ref_df,
  title_text = "Syntenic (genome)",
  x_title = "Reference genome (%)",
  show_y = TRUE,
  xlim_max = SYNTENIC_XMAX
)

p_syn_query_label <- make_metric_plot(
  metric_query_df,
  title_text = "Syntenic (query)",
  x_title = "Query genome (%)",
  show_y = TRUE,
  xlim_max = SYNTENIC_XMAX
)

p_comp_genome_label <- make_comp_plot(
  comp_ref_df,
  title_text = "Composition (genome)",
  x_title = "Reference genome (%)",
  show_y = TRUE,
  xlim_max = COMPOSITION_XMAX
)

p_comp_query_label <- make_comp_plot(
  comp_query_df,
  title_text = "Composition (query)",
  x_title = "Query genome (%)",
  show_y = TRUE,
  xlim_max = COMPOSITION_XMAX
)

# =========================================================
# No-label plots
# y-axis tick marks are still kept:
# one tick per species, but no species names.
# =========================================================
p_syn_genome_nolabel <- make_metric_plot(
  metric_ref_df,
  title_text = "Syntenic (genome)",
  x_title = "Reference genome (%)",
  show_y = FALSE,
  xlim_max = SYNTENIC_XMAX
)

p_syn_query_nolabel <- make_metric_plot(
  metric_query_df,
  title_text = "Syntenic (query)",
  x_title = "Query genome (%)",
  show_y = FALSE,
  xlim_max = SYNTENIC_XMAX
)

p_comp_genome_nolabel <- make_comp_plot(
  comp_ref_df,
  title_text = "Composition (genome)",
  x_title = "Reference genome (%)",
  show_y = FALSE,
  xlim_max = COMPOSITION_XMAX
)

p_comp_query_nolabel <- make_comp_plot(
  comp_query_df,
  title_text = "Composition (query)",
  x_title = "Query genome (%)",
  show_y = FALSE,
  xlim_max = COMPOSITION_XMAX
)

# =========================================================
# Combined plots
# =========================================================
p_syn_combo_label <- (p_syn_genome_label | p_syn_query_label) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

p_syn_combo_nolabel <- (p_syn_genome_nolabel | p_syn_query_nolabel) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

p_comp_combo_label <- (p_comp_genome_label | p_comp_query_label) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

p_comp_combo_nolabel <- (p_comp_genome_nolabel | p_comp_query_nolabel) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

# =========================================================
# Save
# =========================================================
save_plot <- function(p, file_stub, width, height) {
  ggsave(
    filename = file.path(outdir, paste0(file_stub, ".pdf")),
    plot = p,
    width = width,
    height = height,
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(outdir, paste0(file_stub, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )
}

save_plot(p_syn_genome_label, "01_syntenic_genome_grouped_label", WIDTH_LABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_syn_query_label,  "02_syntenic_query_grouped_label",  WIDTH_LABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_syn_combo_label,  "03_syntenic_grouped_2panel_label", WIDTH_LABEL_COMBO, HEIGHT_COMBO)

save_plot(p_comp_genome_label, "04_composition_genome_grouped_label", WIDTH_LABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_comp_query_label,  "05_composition_query_grouped_label",  WIDTH_LABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_comp_combo_label,  "06_composition_grouped_2panel_label", WIDTH_LABEL_COMBO, HEIGHT_COMBO)

save_plot(p_syn_genome_nolabel, "07_syntenic_genome_grouped_nolabel", WIDTH_NOLABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_syn_query_nolabel,  "08_syntenic_query_grouped_nolabel",  WIDTH_NOLABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_syn_combo_nolabel,  "09_syntenic_grouped_2panel_nolabel", WIDTH_NOLABEL_COMBO, HEIGHT_COMBO)

save_plot(p_comp_genome_nolabel, "10_composition_genome_grouped_nolabel", WIDTH_NOLABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_comp_query_nolabel,  "11_composition_query_grouped_nolabel",  WIDTH_NOLABEL_SINGLE, HEIGHT_SINGLE)
save_plot(p_comp_combo_nolabel,  "12_composition_grouped_2panel_nolabel", WIDTH_NOLABEL_COMBO, HEIGHT_COMBO)

message("Done. Output directory: ", outdir)
message("Species slots from tree/order file kept: ", nrow(species_meta))
message("Each species has three reference slots: ", paste(target_refs_raw, collapse = ", "))
message("Empty species in tree, such as Coix_aquatica if no OK result exists, are retained as blank vertical space.")
message("Y-axis ticks: one tick per species, not one tick per bar.")
message("Self-reference values removed from bars: ", REMOVE_SELF_VALUES)
message("Syntenic x-axis upper limit fixed at: 95%")
message("Grid: light gray solid vertical lines only")
