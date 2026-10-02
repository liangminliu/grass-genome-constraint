#!/usr/bin/env Rscript

## =========================================================
## Unified gene-class SNP density, deleterious burden density,
## deleterious ratio, and CDS GERP constrained-site percentage
## across moso bamboo, teosinte/maize, and Tausch's goatgrass.
##
## Main design changes in this version:
##   1) SNP density is drawn as a standalone species-combined row:
##      three species panels, one shared Y-axis title, free Y ticks.
##   2) All panels use the same compact layout and numerical x positions,
##      so class groups are closer and panel sizes are smaller.
##   3) Fill colors preserve species color identity:
##      moso bamboo = blue series, teosinte/maize = brown-orange series,
##      Tausch's goatgrass = green series.
##      Within each species, Core/Syntenic/Nonsyntenic are represented by
##      dark/middle/light shades. The figure legend uses grey shades only
##      to explain the class shade order.
##   4) Long violin tails are handled by display-only q99 clipping.
##      All statistics use the full data before display clipping.
##   5) Optional ggbreak axis break can be enabled if the package is
##      installed, but q99 clipping is the default and most stable mode.
##
## Required R packages:
##   dplyr, tidyr, ggplot2, patchwork, grid
## Optional:
##   ggbreak, if USE_GGBREAK <- TRUE
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

## =========================================================
## 1. User parameters
## =========================================================

PAN_FILE_CANDIDATES <- c(
  "core_5sp_1to1.SG.pan",
  "../../../01.OG_gene/core_5sp_1to1.SG.pan"
)

THRESHOLDS <- c("gt2", "gt4", "gt6")
DATA_TYPES <- c("ALL", "SNPpos", "SNP0")
PLOT_DATA_TYPE <- "ALL"
SNP_DENSITY_THRESHOLD <- "gt4"

OUTDIR    <- "gene_class_unified_SNP_GERP_A4_all_modes_results"
TABLE_DIR <- file.path(OUTDIR, "tables")
STAT_DIR  <- file.path(OUTDIR, "stats")
FIG_DIR   <- file.path(OUTDIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

## Plot geometry mode: "box" or "violin_box".
## The script writes both modes by default.
GEOM_MODES_TO_SAVE <- c("box", "violin_box")

## Display-only tail handling.
DISPLAY_UPPER_QUANTILE <- 0.99
DISPLAY_LOWER_QUANTILE_LOG10 <- 0.01
PSEUDO <- 1e-6

## Optional true broken-axis mode. Default FALSE because q99 clipping is
## safer with significance brackets and does not require extra packages.
USE_GGBREAK <- FALSE
GGBREAK_MIN_TAIL_RATIO <- 2.0

## Significance display.
show_ns_lines <- FALSE
SHOW_SIGNIF <- TRUE
SHOW_N_IN_BOXPLOT <- FALSE

## Compact figure dimensions.
ROW_FIG_WIDTH  <- 7.60
ROW_FIG_HEIGHT <- 2.30
A4_WIDTH       <- 8.27
A4_HEIGHT      <- 11.69
SINGLE_SPECIES_WIDTH  <- 8.00
SINGLE_SPECIES_HEIGHT <- 10.30
FIG_DPI        <- 600

## Class order.
class_order <- c("Core", "Syntenic", "Nonsyntenic")

class_full_names <- c(
  "Core"        = "strict_perfect_copy_syntenic_core",
  "Syntenic"    = "broad_syntenic_excluding_strict_core",
  "Nonsyntenic" = "candidate_nonsyntenic_private"
)

## Species order and labels.
species_order <- c("ped", "Zmay", "Atau")
species_display <- c(
  "ped"  = "Moso bamboo",
  "Zmay" = "Teosinte",
  "Atau" = "Tausch's goatgrass"
)

species_axis_display <- c(
  "ped"  = "italic(P.~edulis)",
  "Zmay" = "italic(Z.~mays)",
  "Atau" = "italic(A.~tauschii)"
)

## Species-based class palettes.
## Core/Syntenic/Nonsyntenic = dark/middle/light of each species color.
species_class_fill <- list(
  ped = c(
    "Core"        = "#2F74B8",
    "Syntenic"    = "#5A97D0",
    "Nonsyntenic" = "#A0C8E8"
  ),
  Zmay = c(
    "Core"        = "#5F3B2E",
    "Syntenic"    = "#8F5D44",
    "Nonsyntenic" = "#CFB8A9"
  ),
  Atau = c(
    "Core"        = "#507D39",
    "Syntenic"    = "#8CB26C",
    "Nonsyntenic" = "#D7E9BC"
  )
)

species_class_line <- list(
  ped = c(
    "Core"        = "#174D83",
    "Syntenic"    = "#2F74B8",
    "Nonsyntenic" = "#5A97D0"
  ),
  Zmay = c(
    "Core"        = "#5F3B2E",
    "Syntenic"    = "#8F5D44",
    "Nonsyntenic" = "#CFB8A9"
  ),
  Atau = c(
    "Core"        = "#2F4E22",
    "Syntenic"    = "#507D39",
    "Nonsyntenic" = "#8CB26C"
  )
)

## Grey legend for gene-class shade interpretation.
class_legend_fill <- c(
  "Core"        = "#4D4D4D",
  "Syntenic"    = "#9A9A9A",
  "Nonsyntenic" = "#D4D4D4"
)

## Numerical offsets. Smaller values make class groups tighter.
## Numerical positions. Smaller offsets and compressed species spacing
## make species groups closer in the A4 combined panels.
## Geometry model:
##   - Boxplot class spacing is controlled by BOX_CLASS_STEP.
##   - Violin class spacing is exactly one third of the boxplot class spacing.
##   - Between two neighboring species groups, the visual gap is one geometry width.
## These values are used dynamically according to geom_mode.
BOX_WIDTH <- 0.078
BOX_CLASS_STEP <- 0.108
VIOLIN_CLASS_STEP <- BOX_CLASS_STEP / 3
VIOLIN_WIDTH <- VIOLIN_CLASS_STEP * 0.78
VIOLIN_BOX_WIDTH <- VIOLIN_WIDTH * 0.40
SPECIES_GAP_WIDTH_MULT <- 1.00

## Kept only for backward compatibility. Current code uses get_species_centers().
SPECIES_CENTER_SPACING <- 2 * BOX_CLASS_STEP + (1 + SPECIES_GAP_WIDTH_MULT) * BOX_WIDTH
class_offsets <- c(
  "Core"        = -BOX_CLASS_STEP,
  "Syntenic"    =  0.000,
  "Nonsyntenic" =  BOX_CLASS_STEP
)

## Vertical layout of the species-local (combined) panels, in normalized
## 0..1 display coordinates. The data band is vertically CENTERED so the
## boxes/violins sit in the middle of each panel. The three vertical spines
## (left border = P. edulis, first separator = Z. mays, right border =
## A. tauschii) are drawn across the FULL panel height as real Y axes, with
## ticks and numbers pointing outward.
LOCAL_BAND_BOTTOM <- 0.170
LOCAL_BAND_TOP    <- 0.680
LOCAL_SIG_BASE    <- 0.740
LOCAL_SIG_STEP    <- 0.048
LOCAL_N_Y         <- 0.070
## Bottom margin brings the spine floor down to the axis baseline (band
## bottom minus this, clamped at 0). Top margin is the small headroom kept
## above the highest significance bracket.
LOCAL_BOTTOM_MARGIN <- 0.170
LOCAL_TOP_MARGIN    <- 0.045
LOCAL_AXIS_OFFSET <- 0.23
LOCAL_AXIS_TICK   <- 0.014
N_LOCAL_AXIS_BREAKS <- 4
FIRST_SEPARATOR_BUFFER <- 0.07

## =========================================================
## 2. Basic helpers
## =========================================================

find_first_existing <- function(paths, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop("Missing ", label, ". Tried: ", paste(paths, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

PAN_FILE <- find_first_existing(PAN_FILE_CANDIDATES, "core_5sp_1to1.SG.pan")

read_tsv_base <- function(fp, header = TRUE) {
  if (!file.exists(fp)) {
    stop("Missing input file: ", fp, call. = FALSE)
  }
  read.delim(
    fp,
    header = header,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv_base <- function(df, fp) {
  write.table(
    df,
    file = fp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

safe_ggsave <- function(filename, plot, width, height, dpi = NULL) {
  if (grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    ok <- FALSE
    try({
      ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = cairo_pdf,
        limitsize = FALSE
      )
      ok <- TRUE
    }, silent = TRUE)
    if (!ok) {
      ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
      )
    }
  } else {
    ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE
    )
  }
}

get_species_label <- function(sp) {
  sp_id <- as.character(sp)[1]
  if (is.na(sp_id) || !(sp_id %in% names(species_display))) {
    stop("Unknown species id: ", paste(sp, collapse = ", "), call. = FALSE)
  }
  unname(species_display[sp_id])
}

strip_prefix <- function(x, sp) {
  sp_id <- as.character(sp)[1]
  x <- as.character(x)

  if (sp_id == "ped") {
    x <- sub("^pedC_", "", x)
    x <- sub("^pedD_", "", x)
    x <- sub("^P\\.edulis_", "", x)
  } else if (sp_id == "Zmay") {
    x <- sub("^Zmay_", "", x)
    x <- sub("^Zea_mays_", "", x)
  } else if (sp_id == "Atau") {
    x <- sub("^Atau_", "", x)
    x <- sub("^Aegilops_tauschii_", "", x)
  } else {
    stop("Unknown species in strip_prefix(): ", sp_id, call. = FALSE)
  }

  x
}

p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

nice_y_upper <- function(x) {
  if (!is.finite(x) || x <= 0) return(1)
  if (x <= 0.001) return(signif(x, 2))
  if (x <= 0.01)  return(ceiling(x / 0.001) * 0.001)
  if (x <= 0.05)  return(ceiling(x / 0.005) * 0.005)
  if (x <= 0.10)  return(ceiling(x / 0.01) * 0.01)
  if (x <= 1.00)  return(ceiling(x / 0.10) * 0.10)
  if (x <= 5)     return(ceiling(x))
  if (x <= 20)    return(ceiling(x / 2) * 2)
  if (x <= 50)    return(ceiling(x / 5) * 5)
  ceiling(x / 10) * 10
}

make_fill_key <- function(sp, cls) {
  paste(as.character(sp), as.character(cls), sep = "__")
}

all_fill_values <- unlist(lapply(names(species_class_fill), function(sp) {
  vals <- species_class_fill[[sp]]
  names(vals) <- make_fill_key(sp, names(vals))
  vals
}), use.names = TRUE)

all_line_values <- unlist(lapply(names(species_class_line), function(sp) {
  vals <- species_class_line[[sp]]
  names(vals) <- make_fill_key(sp, names(vals))
  vals
}), use.names = TRUE)

## =========================================================
## 3. Input file builders
## =========================================================

snp_core_file <- function(sp, thr) {
  sp_id <- as.character(sp)[1]
  if (sp_id == "ped") {
    stop("Use pedC/pedD core builders for ped.", call. = FALSE)
  }
  paste0(sp_id, "_", thr, "_core_gene_deleterious.tsv")
}

snp_pedC_core_file <- function(thr) paste0("pedC_", thr, "_core_gene_deleterious.tsv")
snp_pedD_core_file <- function(thr) paste0("pedD_", thr, "_core_gene_deleterious.tsv")

combined_file <- function(sp, thr) {
  sp_id <- as.character(sp)[1]
  if (sp_id == "ped") return(paste0("ped.gerp_", thr, "_combined.tsv"))
  paste0(sp_id, "_gerp_", thr, "_combined.tsv")
}

gerp_core_per_gene_file <- function(sp, thr) {
  sp_id <- as.character(sp)[1]
  if (sp_id == "ped") {
    stop("Use pedC/pedD GERP per-gene builders for ped.", call. = FALSE)
  }
  paste0(sp_id, "_gerp_", thr, "_CDS_per_gene.tsv")
}

gerp_pedC_file <- function(thr) paste0("pedC_gerp_", thr, "_CDS_per_gene.tsv")
gerp_pedD_file <- function(thr) paste0("pedD_gerp_", thr, "_CDS_per_gene.tsv")

other_file <- function(sp) paste0(as.character(sp)[1], "_other_syntenic.genes.txt")
private_file <- function(sp) paste0(as.character(sp)[1], "_candidate_nonsyntenic_private.genes.txt")

## =========================================================
## 4. Gene-class maps
## =========================================================

read_gene_list <- function(fp, sp) {
  sp_id <- as.character(sp)[1]
  if (!file.exists(fp)) {
    stop("Missing gene list file: ", fp, call. = FALSE)
  }
  x <- scan(fp, what = character(), quiet = TRUE)
  x <- strip_prefix(x, sp_id)
  x <- x[!is.na(x) & x != "" & x != "-"]
  unique(x)
}

read_pan <- function() {
  pan <- read_tsv_base(PAN_FILE, header = FALSE)
  if (ncol(pan) < 9) {
    stop("core_5sp_1to1.SG.pan should contain at least 9 columns.", call. = FALSE)
  }
  pan
}

read_core_genes_nonped <- function(sp) {
  sp_id <- as.character(sp)[1]
  pan <- read_pan()
  col_id <- if (sp_id == "Atau") 5 else if (sp_id == "Zmay") 9 else {
    stop("read_core_genes_nonped() only supports Atau and Zmay.", call. = FALSE)
  }
  genes <- strip_prefix(pan[[col_id]], sp_id)
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != "" & genes != "-"]
  cat("[OK] Strict core ", sp_id, " genes from SG pan: ", length(genes), "\n", sep = "")
  genes
}

build_gene_class_map_nonped <- function(sp) {
  sp_id <- as.character(sp)[1]
  core_genes <- read_core_genes_nonped(sp_id)
  other_genes <- read_gene_list(other_file(sp_id), sp_id)
  private_genes <- read_gene_list(private_file(sp_id), sp_id)

  other_clean <- setdiff(other_genes, core_genes)
  private_clean <- setdiff(private_genes, unique(c(core_genes, other_clean)))

  out <- bind_rows(
    data.frame(unit_id = core_genes, gene_class = "Core", stringsAsFactors = FALSE),
    data.frame(unit_id = other_clean, gene_class = "Syntenic", stringsAsFactors = FALSE),
    data.frame(unit_id = private_clean, gene_class = "Nonsyntenic", stringsAsFactors = FALSE)
  ) %>%
    distinct(unit_id, .keep_all = TRUE) %>%
    mutate(
      species = sp_id,
      species_label = get_species_label(sp_id),
      gene_class = factor(gene_class, levels = class_order)
    )

  cat("[OK] ", sp_id, " gene-class map:\n", sep = "")
  print(table(out$gene_class, useNA = "ifany"))
  out
}

read_core_sg_map_ped <- function() {
  pan <- read_pan()
  out <- pan[, 1:9]
  colnames(out) <- c("SG", "n_species", "n_genes", "flag", "Atau", "Osat", "pedC", "pedD", "Zmay")

  out %>%
    transmute(
      SG = as.character(SG),
      pedC = strip_prefix(pedC, "ped"),
      pedD = strip_prefix(pedD, "ped")
    ) %>%
    filter(
      !is.na(SG),
      !is.na(pedC), pedC != "", pedC != "-",
      !is.na(pedD), pedD != "", pedD != "-"
    ) %>%
    distinct(SG, .keep_all = TRUE)
}

build_gene_class_map_ped_genes_for_gerp <- function() {
  pan <- read_pan()
  core_genes <- unique(c(strip_prefix(pan[[7]], "ped"), strip_prefix(pan[[8]], "ped")))
  core_genes <- core_genes[!is.na(core_genes) & core_genes != "" & core_genes != "-"]

  other_genes <- read_gene_list(other_file("ped"), "ped")
  private_genes <- read_gene_list(private_file("ped"), "ped")

  other_clean <- setdiff(other_genes, core_genes)
  private_clean <- setdiff(private_genes, unique(c(core_genes, other_clean)))

  out <- bind_rows(
    data.frame(unit_id = core_genes, gene_class = "Core", stringsAsFactors = FALSE),
    data.frame(unit_id = other_clean, gene_class = "Syntenic", stringsAsFactors = FALSE),
    data.frame(unit_id = private_clean, gene_class = "Nonsyntenic", stringsAsFactors = FALSE)
  ) %>%
    distinct(unit_id, .keep_all = TRUE) %>%
    mutate(
      species = "ped",
      species_label = get_species_label("ped"),
      gene_class = factor(gene_class, levels = class_order)
    )

  cat("[OK] ped GERP gene-class map:\n")
  print(table(out$gene_class, useNA = "ifany"))
  out
}

## =========================================================
## 5. SNP/deleterious metrics
## =========================================================

read_snp_core_stats_nonped <- function(sp, thr, map_obj) {
  sp_id <- as.character(sp)[1]
  fp <- snp_core_file(sp_id, thr)
  df <- read_tsv_base(fp, header = TRUE)

  required_cols <- c("gene_id", "cds_len", "cds_snp_total", "cds_deleterious_snp")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }

  core_units <- map_obj %>% filter(gene_class == "Core") %>% pull(unit_id)

  df %>%
    transmute(
      species = sp_id,
      species_label = get_species_label(sp_id),
      threshold = thr,
      gene_class = "Core",
      unit_id = strip_prefix(.data$gene_id, sp_id),
      cds_len = suppressWarnings(as.numeric(.data$cds_len)),
      cds_snp_total = suppressWarnings(as.numeric(.data$cds_snp_total)),
      cds_deleterious_snp = suppressWarnings(as.numeric(.data$cds_deleterious_snp)),
      source = fp
    ) %>%
    filter(
      unit_id %in% core_units,
      !is.na(cds_len), cds_len > 0,
      !is.na(cds_snp_total),
      !is.na(cds_deleterious_snp)
    )
}

read_snp_noncore_stats_nonped <- function(sp, thr, map_obj) {
  sp_id <- as.character(sp)[1]
  fp <- combined_file(sp_id, thr)
  df <- read_tsv_base(fp, header = TRUE)

  required_cols <- c("group", "gene_id", "cds_len", "cds_snp_total", "cds_deleterious_snp")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }

  noncore_map <- map_obj %>% filter(gene_class %in% c("Syntenic", "Nonsyntenic"))

  df %>%
    mutate(
      unit_id = strip_prefix(.data$gene_id, sp_id),
      group_class = case_when(
        .data$group == "other_syntenic" ~ "Syntenic",
        .data$group == "candidate_nonsyntenic_private" ~ "Nonsyntenic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group_class)) %>%
    inner_join(noncore_map %>% select(unit_id, gene_class), by = "unit_id") %>%
    filter(as.character(.data$gene_class) == .data$group_class) %>%
    transmute(
      species = sp_id,
      species_label = get_species_label(sp_id),
      threshold = thr,
      gene_class = as.character(.data$gene_class),
      unit_id = .data$unit_id,
      cds_len = suppressWarnings(as.numeric(.data$cds_len)),
      cds_snp_total = suppressWarnings(as.numeric(.data$cds_snp_total)),
      cds_deleterious_snp = suppressWarnings(as.numeric(.data$cds_deleterious_snp)),
      source = paste0(fp, ":", .data$group)
    ) %>%
    filter(
      !is.na(cds_len), cds_len > 0,
      !is.na(cds_snp_total),
      !is.na(cds_deleterious_snp)
    )
}

read_snp_core_subgenome_stats_ped <- function(fp, label) {
  df <- read_tsv_base(fp, header = TRUE)
  required_cols <- c("gene_id", "cds_len", "cds_snp_total", "cds_deleterious_snp")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df %>%
    transmute(
      gene_id = strip_prefix(.data$gene_id, "ped"),
      cds_len = suppressWarnings(as.numeric(.data$cds_len)),
      cds_snp_total = suppressWarnings(as.numeric(.data$cds_snp_total)),
      cds_deleterious_snp = suppressWarnings(as.numeric(.data$cds_deleterious_snp))
    ) %>%
    filter(
      !is.na(gene_id), gene_id != "",
      !is.na(cds_len), cds_len > 0,
      !is.na(cds_snp_total),
      !is.na(cds_deleterious_snp)
    ) %>%
    group_by(gene_id) %>%
    summarise(
      cds_len = mean(cds_len, na.rm = TRUE),
      cds_snp_total = mean(cds_snp_total, na.rm = TRUE),
      cds_deleterious_snp = mean(cds_deleterious_snp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(label, "_", .x), c(cds_len, cds_snp_total, cds_deleterious_snp))
}

read_snp_ped_core_sg <- function(thr, core_sg_map) {
  pedC_df <- read_snp_core_subgenome_stats_ped(snp_pedC_core_file(thr), "pedC")
  pedD_df <- read_snp_core_subgenome_stats_ped(snp_pedD_core_file(thr), "pedD")

  core_sg_map %>%
    left_join(pedC_df, by = c("pedC" = "gene_id")) %>%
    left_join(pedD_df, by = c("pedD" = "gene_id")) %>%
    filter(!is.na(.data$pedC_cds_len), !is.na(.data$pedD_cds_len)) %>%
    transmute(
      species = "ped",
      species_label = get_species_label("ped"),
      threshold = thr,
      gene_class = "Core",
      unit_id = .data$SG,
      pedC_gene = .data$pedC,
      pedD_gene = .data$pedD,
      cds_len = .data$pedC_cds_len + .data$pedD_cds_len,
      cds_snp_total = .data$pedC_cds_snp_total + .data$pedD_cds_snp_total,
      cds_deleterious_snp = .data$pedC_cds_deleterious_snp + .data$pedD_cds_deleterious_snp,
      source = "core_SG_pedC_plus_pedD"
    )
}

read_snp_ped_noncore <- function(thr) {
  fp <- combined_file("ped", thr)
  df <- read_tsv_base(fp, header = TRUE)
  required_cols <- c("group", "gene_id", "cds_len", "cds_snp_total", "cds_deleterious_snp")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df %>%
    mutate(
      gene_class = case_when(
        .data$group == "other_syntenic" ~ "Syntenic",
        .data$group == "candidate_nonsyntenic_private" ~ "Nonsyntenic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(gene_class)) %>%
    transmute(
      species = "ped",
      species_label = get_species_label("ped"),
      threshold = thr,
      gene_class = .data$gene_class,
      unit_id = strip_prefix(.data$gene_id, "ped"),
      pedC_gene = NA_character_,
      pedD_gene = NA_character_,
      cds_len = suppressWarnings(as.numeric(.data$cds_len)),
      cds_snp_total = suppressWarnings(as.numeric(.data$cds_snp_total)),
      cds_deleterious_snp = suppressWarnings(as.numeric(.data$cds_deleterious_snp)),
      source = paste0("combined_", .data$group)
    ) %>%
    filter(
      !is.na(unit_id), unit_id != "",
      !is.na(cds_len), cds_len > 0,
      !is.na(cds_snp_total),
      !is.na(cds_deleterious_snp)
    )
}

build_snp_one_threshold_species <- function(sp, thr, map_obj) {
  sp_id <- as.character(sp)[1]
  if (sp_id == "ped") {
    core_df <- read_snp_ped_core_sg(thr, map_obj)
    noncore_df <- read_snp_ped_noncore(thr)
  } else {
    core_df <- read_snp_core_stats_nonped(sp_id, thr, map_obj)
    noncore_df <- read_snp_noncore_stats_nonped(sp_id, thr, map_obj)
  }

  out <- bind_rows(core_df, noncore_df) %>%
    mutate(
      species = factor(.data$species, levels = species_order),
      species_label = factor(.data$species_label, levels = unname(species_display[species_order])),
      threshold = factor(.data$threshold, levels = THRESHOLDS),
      gene_class = factor(.data$gene_class, levels = class_order),
      cds_snp_density = .data$cds_snp_total / .data$cds_len,
      cds_deleterious_ratio = ifelse(.data$cds_snp_total > 0, .data$cds_deleterious_snp / .data$cds_snp_total, 0),
      cds_deleterious_density = .data$cds_deleterious_snp / .data$cds_len
    ) %>%
    filter(
      !is.na(.data$gene_class),
      is.finite(.data$cds_snp_density),
      is.finite(.data$cds_deleterious_ratio),
      is.finite(.data$cds_deleterious_density)
    )

  cat("[INFO] SNP metrics: species=", sp_id, ", threshold=", thr, "\n", sep = "")
  print(table(out$gene_class, useNA = "ifany"))
  out
}

split_by_snp_class <- function(df) {
  list(
    ALL = df,
    SNPpos = df %>% filter(.data$cds_snp_total > 0),
    SNP0 = df %>% filter(.data$cds_snp_total == 0)
  )
}

run_stats_metric <- function(df_plot, metric_col, metric_label, data_type = NA_character_) {
  df_plot <- df_plot %>%
    filter(!is.na(.data$species), !is.na(.data$gene_class), is.finite(.data[[metric_col]])) %>%
    droplevels()

  out <- list()
  idx <- 1

  for (sp_id in intersect(species_order, unique(as.character(df_plot$species)))) {
    for (thr in intersect(THRESHOLDS, unique(as.character(df_plot$threshold)))) {
      sub_df <- df_plot %>%
        filter(as.character(.data$species) == sp_id, as.character(.data$threshold) == thr) %>%
        droplevels()

      observed_levels <- levels(droplevels(sub_df$gene_class))

      if (length(observed_levels) < 2 || nrow(sub_df) < 3) {
        out[[idx]] <- tibble(
          species = sp_id,
          species_label = get_species_label(sp_id),
          threshold = thr,
          data_type = data_type,
          metric = metric_label,
          value_col = metric_col,
          test = NA_character_,
          KW_chisq = NA_real_,
          KW_df = NA_real_,
          KW_p = NA_real_,
          Comparison = NA_character_,
          group1 = NA_character_,
          group2 = NA_character_,
          n_group1 = NA_integer_,
          n_group2 = NA_integer_,
          median_group1 = NA_real_,
          median_group2 = NA_real_,
          median_diff_group1_minus_group2 = NA_real_,
          P.raw = NA_real_,
          P.adj = NA_real_,
          label = NA_character_
        )
        idx <- idx + 1
        next
      }

      kw <- kruskal.test(sub_df[[metric_col]] ~ sub_df$gene_class)
      pair_combs <- combn(observed_levels, 2, simplify = FALSE)

      pair_tbl <- lapply(pair_combs, function(cc) {
        g1 <- cc[1]
        g2 <- cc[2]
        x <- sub_df[[metric_col]][sub_df$gene_class == g1]
        y <- sub_df[[metric_col]][sub_df$gene_class == g2]

        wt <- tryCatch(
          wilcox.test(x, y, paired = FALSE, exact = FALSE),
          error = function(e) NULL
        )

        tibble(
          species = sp_id,
          species_label = get_species_label(sp_id),
          threshold = thr,
          data_type = data_type,
          metric = metric_label,
          value_col = metric_col,
          test = "Kruskal-Wallis + Wilcoxon rank-sum; BH correction",
          KW_chisq = unname(kw$statistic),
          KW_df = unname(kw$parameter),
          KW_p = kw$p.value,
          Comparison = paste0(class_full_names[[g1]], " - ", class_full_names[[g2]]),
          group1 = class_full_names[[g1]],
          group2 = class_full_names[[g2]],
          n_group1 = length(x),
          n_group2 = length(y),
          median_group1 = median(x, na.rm = TRUE),
          median_group2 = median(y, na.rm = TRUE),
          median_diff_group1_minus_group2 = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
          P.raw = ifelse(is.null(wt), NA_real_, wt$p.value)
        )
      }) %>%
        bind_rows() %>%
        mutate(
          P.adj = p.adjust(.data$P.raw, method = "BH"),
          label = vapply(.data$P.adj, p_to_star, character(1))
        )

      out[[idx]] <- pair_tbl
      idx <- idx + 1
    }
  }

  bind_rows(out)
}

summarise_metric_table <- function(df, data_type_label = NA_character_) {
  ## Important:
  ## Do not name this argument `data_type`. Some input tables already contain
  ## a `data_type` column. In dplyr data-masked verbs, `data_type` would then be
  ## resolved as the whole column instead of a scalar function argument, causing
  ## errors such as: `data_type must be size 1, not n`.
  dt_label <- as.character(data_type_label)[1]
  if (length(dt_label) == 0 || is.na(dt_label)) dt_label <- NA_character_

  has_cds_snp_total <- "cds_snp_total" %in% names(df)
  has_cds_len <- "cds_len" %in% names(df)
  has_snp_density <- "cds_snp_density" %in% names(df)
  has_del_density <- "cds_deleterious_density" %in% names(df)
  has_del_ratio <- "cds_deleterious_ratio" %in% names(df)
  has_gerp_pct <- "cds_conserved_pct" %in% names(df)

  df %>%
    group_by(.data$species, .data$species_label, .data$threshold, .data$gene_class) %>%
    summarise(
      data_type = .env$dt_label,
      n_units = dplyr::n(),
      n_snp_pos = if (.env$has_cds_snp_total) sum(.data$cds_snp_total > 0, na.rm = TRUE) else NA_integer_,
      n_snp0 = if (.env$has_cds_snp_total) sum(.data$cds_snp_total == 0, na.rm = TRUE) else NA_integer_,
      median_cds_len = if (.env$has_cds_len) median(.data$cds_len, na.rm = TRUE) else NA_real_,
      median_cds_snp_density = if (.env$has_snp_density) median(.data$cds_snp_density, na.rm = TRUE) else NA_real_,
      median_cds_deleterious_density = if (.env$has_del_density) median(.data$cds_deleterious_density, na.rm = TRUE) else NA_real_,
      median_cds_deleterious_ratio = if (.env$has_del_ratio) median(.data$cds_deleterious_ratio, na.rm = TRUE) else NA_real_,
      median_cds_conserved_pct = if (.env$has_gerp_pct) median(.data$cds_conserved_pct, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}

## =========================================================
## 6. GERP constrained CDS percentage metrics
## =========================================================

get_cds_pct_column <- function(df, fp) {
  if ("cds_conserved_pct" %in% colnames(df)) {
    return(suppressWarnings(as.numeric(df$cds_conserved_pct)))
  }
  if ("cds_pct" %in% colnames(df)) {
    return(suppressWarnings(as.numeric(df$cds_pct)))
  }
  if (all(c("cds_conserved_sites", "cds_len") %in% colnames(df))) {
    cds_sites <- suppressWarnings(as.numeric(df$cds_conserved_sites))
    cds_len <- suppressWarnings(as.numeric(df$cds_len))
    return(100 * cds_sites / cds_len)
  }
  if (all(c("cds_sites", "cds_len") %in% colnames(df))) {
    cds_sites <- suppressWarnings(as.numeric(df$cds_sites))
    cds_len <- suppressWarnings(as.numeric(df$cds_len))
    return(100 * cds_sites / cds_len)
  }
  stop(
    "Cannot find CDS constrained percentage columns in file: ", fp,
    "\nNeed one of: cds_conserved_pct, cds_pct, cds_conserved_sites/cds_len, or cds_sites/cds_len.",
    call. = FALSE
  )
}

get_cds_len_column <- function(df) {
  if ("cds_len" %in% colnames(df)) {
    return(suppressWarnings(as.numeric(df$cds_len)))
  }
  rep(NA_real_, nrow(df))
}

read_gerp_combined_stats <- function(sp, thr) {
  sp_id <- as.character(sp)[1]
  fp <- combined_file(sp_id, thr)
  if (!file.exists(fp)) {
    warning("Combined GERP file not found: ", fp)
    return(NULL)
  }

  df <- read_tsv_base(fp, header = TRUE)
  if (!"gene_id" %in% colnames(df)) {
    stop("Missing gene_id column in combined file: ", fp, call. = FALSE)
  }

  out <- data.frame(
    species = sp_id,
    species_label = get_species_label(sp_id),
    threshold = thr,
    unit_id = strip_prefix(df$gene_id, sp_id),
    cds_len = get_cds_len_column(df),
    cds_conserved_pct = get_cds_pct_column(df, fp),
    source = "combined",
    source_rank = 1,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(.data$unit_id), .data$unit_id != "", !is.na(.data$cds_conserved_pct))

  if (any(!is.na(out$cds_len))) {
    out <- out %>% filter(!is.na(.data$cds_len), .data$cds_len > 0)
  }

  out
}

read_gerp_core_per_gene_nonped <- function(sp, thr) {
  sp_id <- as.character(sp)[1]
  fp <- gerp_core_per_gene_file(sp_id, thr)
  if (!file.exists(fp)) return(NULL)

  df <- read_tsv_base(fp, header = TRUE)
  if (!"gene_id" %in% colnames(df)) {
    colnames(df)[1] <- "gene_id"
  }

  out <- data.frame(
    species = sp_id,
    species_label = get_species_label(sp_id),
    threshold = thr,
    unit_id = strip_prefix(df$gene_id, sp_id),
    cds_len = get_cds_len_column(df),
    cds_conserved_pct = get_cds_pct_column(df, fp),
    source = paste0("core_per_gene_", thr),
    source_rank = 2,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(.data$unit_id), .data$unit_id != "", !is.na(.data$cds_conserved_pct))

  if (any(!is.na(out$cds_len))) {
    out <- out %>% filter(!is.na(.data$cds_len), .data$cds_len > 0)
  }

  out
}

read_gerp_ped_subgenome_one <- function(fp, label, thr) {
  if (!file.exists(fp)) return(NULL)
  df <- read_tsv_base(fp, header = TRUE)
  if (!"gene_id" %in% colnames(df)) colnames(df)[1] <- "gene_id"

  out <- data.frame(
    species = "ped",
    species_label = get_species_label("ped"),
    threshold = thr,
    unit_id = strip_prefix(df$gene_id, "ped"),
    cds_len = get_cds_len_column(df),
    cds_conserved_pct = get_cds_pct_column(df, fp),
    source = label,
    source_rank = 2,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(.data$unit_id), .data$unit_id != "", !is.na(.data$cds_conserved_pct))

  if (any(!is.na(out$cds_len))) {
    out <- out %>% filter(!is.na(.data$cds_len), .data$cds_len > 0)
  }

  out
}

build_gerp_one_threshold_species <- function(sp, thr, map_obj) {
  sp_id <- as.character(sp)[1]

  if (sp_id == "ped") {
    all_stats <- bind_rows(
      read_gerp_combined_stats("ped", thr),
      read_gerp_ped_subgenome_one(gerp_pedC_file(thr), paste0("pedC_", thr), thr),
      read_gerp_ped_subgenome_one(gerp_pedD_file(thr), paste0("pedD_", thr), thr)
    )
  } else {
    all_stats <- bind_rows(
      read_gerp_combined_stats(sp_id, thr),
      read_gerp_core_per_gene_nonped(sp_id, thr)
    )
  }

  if (nrow(all_stats) == 0) {
    stop("No usable GERP stats found for species=", sp_id, ", threshold=", thr, call. = FALSE)
  }

  all_stats <- all_stats %>%
    arrange(.data$unit_id, .data$source_rank) %>%
    group_by(.data$unit_id) %>%
    summarise(
      species = first(.data$species),
      species_label = first(.data$species_label),
      threshold = first(.data$threshold),
      cds_len = first(.data$cds_len),
      cds_conserved_pct = first(.data$cds_conserved_pct),
      source = first(.data$source),
      .groups = "drop"
    )

  out <- all_stats %>%
    inner_join(map_obj %>% select(unit_id, gene_class), by = "unit_id") %>%
    mutate(
      species = factor(.data$species, levels = species_order),
      species_label = factor(.data$species_label, levels = unname(species_display[species_order])),
      threshold = factor(.data$threshold, levels = THRESHOLDS),
      gene_class = factor(.data$gene_class, levels = class_order)
    ) %>%
    filter(!is.na(.data$gene_class), is.finite(.data$cds_conserved_pct))

  cat("[INFO] GERP CDS: species=", sp_id, ", threshold=", thr, "\n", sep = "")
  print(table(out$gene_class, useNA = "ifany"))
  out
}


## =========================================================
## 7. Plot helpers
## =========================================================

theme_compact <- function(base_size = 7.4) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.4, face = "plain", hjust = 0.5, margin = margin(b = 2.5)),
      axis.text.x = element_text(size = base_size - 0.35, color = "black", margin = margin(t = 1.6)),
      axis.text.y = element_text(size = base_size - 0.5, color = "black"),
      axis.title.y = element_text(size = base_size + 0.25, color = "black", margin = margin(r = 4.0)),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.28),
      axis.ticks = element_line(color = "black", linewidth = 0.24),
      axis.ticks.length = unit(1.35, "pt"),
      panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.16),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.34),
      legend.position = "none",
      plot.margin = margin(3.5, 8, 3.5, 8)
    )
}

transform_metric_values <- function(x, transform_mode = "raw") {
  if (transform_mode == "raw") return(x)
  if (transform_mode == "log10") return(log10(x + PSEUDO))
  stop("Unknown transform_mode: ", transform_mode, call. = FALSE)
}

get_axis_power <- function(value_col = NULL, transform_mode = "raw") {
  if (!identical(transform_mode, "raw")) return(0)
  if (!is.null(value_col) && value_col %in% c("cds_snp_density", "cds_deleterious_density")) {
    return(-3)
  }
  0
}

make_y_label <- function(base_label, value_col = NULL, transform_mode = "raw") {
  axis_power <- get_axis_power(value_col, transform_mode)
  if (identical(transform_mode, "log10")) {
    return(bquote(log[10](.(base_label) + 10^-6)))
  }
  if (is.finite(axis_power) && axis_power != 0) {
    return(bquote(.(base_label)~"(×"~10^.(axis_power)~")"))
  }
  base_label
}

trim_zeros <- function(x) {
  x <- sub("\\.0+$", "", x)
  x <- sub("(\\.[0-9]*?)0+$", "\\1", x)
  x
}

fmt_axis_num <- function(x, axis_power = 0) {
  x <- as.numeric(x) / (10 ^ axis_power)
  vapply(x, function(v) {
    if (!is.finite(v)) return("")
    if (abs(v) < 1e-12) v <- 0
    av <- abs(v)
    if (abs(v - round(v)) < 1e-8) return(format(round(v), trim = TRUE, scientific = FALSE))
    if (av >= 100) return(trim_zeros(format(round(v, 0), trim = TRUE, scientific = FALSE)))
    if (av >= 10)  return(trim_zeros(format(round(v, 1), trim = TRUE, scientific = FALSE, nsmall = 1)))
    if (av >= 1)   return(trim_zeros(format(round(v, 2), trim = TRUE, scientific = FALSE, nsmall = 2)))
    if (av >= 0.1) return(trim_zeros(format(round(v, 2), trim = TRUE, scientific = FALSE, nsmall = 2)))
    if (av >= 0.01) return(trim_zeros(format(round(v, 3), trim = TRUE, scientific = FALSE, nsmall = 3)))
    trim_zeros(format(signif(v, 2), trim = TRUE, scientific = FALSE))
  }, character(1))
}

nice_axis_breaks <- function(lims, axis_power = 0, n = N_LOCAL_AXIS_BREAKS, include_zero = TRUE) {
  vals <- as.numeric(lims)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(numeric(0))
  lo <- min(vals, na.rm = TRUE)
  hi <- max(vals, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi)) return(numeric(0))
  if (hi <= lo) {
    pad <- max(abs(hi) * 0.05, 1e-8)
    lo <- lo - pad
    hi <- hi + pad
  }

  sf <- 10 ^ axis_power
  lo_s <- lo / sf
  hi_s <- hi / sf
  if (include_zero && lo_s >= 0) lo_s <- 0

  br_s <- pretty(c(lo_s, hi_s), n = n)
  br_s <- br_s[br_s >= lo_s - 1e-9 & br_s <= hi_s + 1e-9]
  if (length(br_s) < 2) br_s <- seq(lo_s, hi_s, length.out = max(2, n))

  unique(br_s * sf)
}

get_display_window <- function(values, transform_mode = "raw", ratio_metric = FALSE,
                               negative_pad = FALSE) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(list(lower = 0, upper = 1, range = 1, base = 1.05, step = 0.08, max_full = 1, upper_core = 1))
  }

  max_full <- max(values, na.rm = TRUE)

  if (transform_mode == "raw") {
    min_val <- min(values, na.rm = TRUE)
    upper_core <- as.numeric(quantile(values, probs = DISPLAY_UPPER_QUANTILE, na.rm = TRUE, names = FALSE))
    lower <- ifelse(min_val >= 0, 0, min_val)

    if (negative_pad && min_val >= 0) {
      pad_mag <- max(upper_core * 0.45, max_full * 0.08, 0.002)
      lower <- -pad_mag
    }
    if (ratio_metric) upper_core <- min(max(upper_core, 0.05), 1)
  } else {
    lower <- as.numeric(quantile(values, probs = DISPLAY_LOWER_QUANTILE_LOG10, na.rm = TRUE, names = FALSE))
    upper_core <- as.numeric(quantile(values, probs = DISPLAY_UPPER_QUANTILE, na.rm = TRUE, names = FALSE))
  }

  if (!is.finite(lower)) lower <- min(values, na.rm = TRUE)
  if (!is.finite(upper_core)) upper_core <- max(values, na.rm = TRUE)
  if (!is.finite(lower)) lower <- 0
  if (!is.finite(upper_core)) upper_core <- 1

  if (upper_core <= lower) {
    pad <- max(abs(upper_core) * 0.10, ifelse(transform_mode == "raw", 0.001, 0.25))
    lower <- upper_core - pad
    upper_core <- upper_core + pad
    if (transform_mode == "raw" && lower < 0 && !negative_pad) lower <- 0
  }

  yrange <- upper_core - lower
  if (!is.finite(yrange) || yrange <= 0) yrange <- max(abs(upper_core), 1)

  step <- max(0.045 * yrange, ifelse(transform_mode == "raw", 0.0003, 0.08))
  base <- upper_core + 0.040 * yrange
  upper <- upper_core + 0.22 * yrange

  list(
    lower = lower,
    upper = upper,
    range = yrange,
    base = base,
    step = step,
    max_full = max_full,
    upper_core = upper_core
  )
}


get_geom_config <- function(geom_mode = "box") {
  if (geom_mode == "box") {
    class_step <- BOX_CLASS_STEP
    geom_width <- BOX_WIDTH
    inner_width <- BOX_WIDTH
  } else if (geom_mode == "violin_box") {
    class_step <- VIOLIN_CLASS_STEP
    geom_width <- VIOLIN_WIDTH
    inner_width <- VIOLIN_BOX_WIDTH
  } else {
    stop("Unknown geom_mode: ", geom_mode, call. = FALSE)
  }

  species_center_spacing <- 2 * class_step + (1 + SPECIES_GAP_WIDTH_MULT) * geom_width
  class_offsets_local <- c(
    "Core"        = -class_step,
    "Syntenic"    =  0,
    "Nonsyntenic" =  class_step
  )

  list(
    class_step = class_step,
    geom_width = geom_width,
    inner_width = inner_width,
    species_center_spacing = species_center_spacing,
    class_offsets = class_offsets_local
  )
}

get_species_centers <- function(species_in_panel, geom_mode = "box") {
  species_in_panel <- as.character(species_in_panel)
  cfg <- get_geom_config(geom_mode)
  center_values <- 1 + (seq_along(species_in_panel) - 1) * cfg$species_center_spacing
  if (length(center_values) >= 2) {
    center_values[1] <- center_values[1] - FIRST_SEPARATOR_BUFFER / 2
    center_values[2] <- center_values[2] + FIRST_SEPARATOR_BUFFER / 2
  }
  setNames(center_values, species_in_panel)
}

add_numeric_positions <- function(df, species_in_panel, geom_mode = "box") {
  species_in_panel <- as.character(species_in_panel)
  cfg <- get_geom_config(geom_mode)
  center_map <- get_species_centers(species_in_panel, geom_mode)

  df %>%
    mutate(
      species_chr = as.character(.data$species),
      gene_class_chr = as.character(.data$gene_class),
      x_center = unname(center_map[.data$species_chr]),
      x_offset = unname(cfg$class_offsets[.data$gene_class_chr]),
      xpos = .data$x_center + .data$x_offset,
      fill_key = make_fill_key(.data$species_chr, .data$gene_class_chr)
    )
}

make_sig_numeric <- function(stat_tbl, species_in_panel, y_base, y_step,
                             local_scale = FALSE, geom_mode = "box") {
  if (!SHOW_SIGNIF) return(tibble())

  sig_df <- stat_tbl %>%
    filter(!is.na(.data$group1), !is.na(.data$group2), !is.na(.data$label))

  if (!show_ns_lines) sig_df <- sig_df %>% filter(.data$label != "ns")
  if (nrow(sig_df) == 0) return(sig_df[0, ])

  inv_map <- setNames(names(class_full_names), class_full_names)
  cfg <- get_geom_config(geom_mode)
  center_map <- get_species_centers(species_in_panel, geom_mode)

  sig_df <- sig_df %>%
    mutate(
      species_chr = as.character(.data$species),
      group1_plot = inv_map[.data$group1],
      group2_plot = inv_map[.data$group2],
      pair_id = paste(pmin(.data$group1_plot, .data$group2_plot), pmax(.data$group1_plot, .data$group2_plot), sep = "__"),
      x_center = unname(center_map[.data$species_chr]),
      x1 = .data$x_center + unname(cfg$class_offsets[.data$group1_plot]),
      x2 = .data$x_center + unname(cfg$class_offsets[.data$group2_plot]),
      xmin = pmin(.data$x1, .data$x2),
      xmax = pmax(.data$x1, .data$x2)
    ) %>%
    group_by(.data$species_chr, .data$pair_id) %>%
    slice_min(order_by = .data$P.adj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(.data$species_chr) %>%
    mutate(
      pair_order = match(.data$pair_id,
                         c("Core__Syntenic", "Core__Nonsyntenic", "Nonsyntenic__Syntenic"))
    ) %>%
    arrange(.data$pair_order, .by_group = TRUE) %>%
    slice_head(n = 3) %>%
    mutate(sig_rank = row_number()) %>%
    ungroup()

  if (local_scale) {
    sig_df <- sig_df %>%
      mutate(
        y = LOCAL_SIG_BASE + (.data$sig_rank - 1) * LOCAL_SIG_STEP,
        tip = 0.012
      )
  } else {
    sig_df <- sig_df %>%
      mutate(
        y = y_base + (.data$sig_rank - 1) * y_step,
        tip = 0.20 * y_step
      )
  }

  sig_df
}

make_n_labels <- function(plot_df, y_lower, y_range, local_scale = FALSE) {
  out <- plot_df %>%
    group_by(.data$species, .data$species_label, .data$gene_class, .data$xpos) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0("n=", .data$n))

  if (local_scale) {
    out <- out %>% mutate(y = LOCAL_N_Y)
  } else {
    out <- out %>% mutate(y = y_lower - 0.090 * y_range)
  }
  out
}

build_local_axis_df <- function(plot_df, species_in_panel, transform_mode = "raw", ratio_metric = FALSE,
                                geom_mode = "box", negative_pad = FALSE, axis_power = 0) {
  win_list <- lapply(species_in_panel, function(sp_id) {
    vals <- plot_df$plot_value[as.character(plot_df$species) == sp_id]
    win <- get_display_window(vals, transform_mode = transform_mode, ratio_metric = ratio_metric,
                              negative_pad = negative_pad)
    data.frame(
      species_chr = sp_id,
      lower = win$lower,
      upper_core = win$upper_core,
      range = win$range,
      center = unname(get_species_centers(species_in_panel, geom_mode)[sp_id]),
      stringsAsFactors = FALSE
    )
  })
  win_df <- bind_rows(win_list)

  axis_df <- bind_rows(lapply(seq_len(nrow(win_df)), function(i) {
    row <- win_df[i, ]
    if (!is.finite(row$lower)) row$lower <- 0
    if (!is.finite(row$upper_core)) row$upper_core <- row$lower + 1
    if (row$upper_core <= row$lower) {
      br <- rep(row$upper_core, N_LOCAL_AXIS_BREAKS)
      y_scaled <- rep((LOCAL_BAND_BOTTOM + LOCAL_BAND_TOP) / 2, N_LOCAL_AXIS_BREAKS)
    } else {
      br <- nice_axis_breaks(
        c(row$lower, row$upper_core),
        axis_power = axis_power,
        n = N_LOCAL_AXIS_BREAKS,
        include_zero = identical(transform_mode, "raw") && row$lower >= 0
      )
      br <- br[br >= row$lower - 1e-12 & br <= row$upper_core + 1e-12]
      if (length(br) < 2) br <- seq(row$lower, row$upper_core, length.out = N_LOCAL_AXIS_BREAKS)
      y_scaled <- LOCAL_BAND_BOTTOM +
        (br - row$lower) / (row$upper_core - row$lower) * (LOCAL_BAND_TOP - LOCAL_BAND_BOTTOM)
    }
    data.frame(
      species_chr = row$species_chr,
      center = row$center,
      tick_value = br,
      y_scaled = y_scaled,
      tick_label = fmt_axis_num(br, axis_power = axis_power),
      stringsAsFactors = FALSE
    )
  }))

  list(win_df = win_df, axis_df = axis_df)
}

plot_numeric_panel <- function(df, stat_tbl, value_col, title, ylab = NULL,
                               species_in_panel = species_order,
                               transform_mode = "raw",
                               geom_mode = "box",
                               show_x_species_labels = TRUE,
                               show_n = FALSE,
                               base_size = 7.4,
                               ratio_metric = FALSE,
                               scale_mode = c("global", "species_local")) {
  scale_mode <- match.arg(scale_mode)
  species_in_panel <- as.character(species_in_panel)
  cfg <- get_geom_config(geom_mode)
  axis_power <- get_axis_power(value_col, transform_mode)
  use_negative_pad <- FALSE

  plot_df <- df %>%
    filter(as.character(.data$species) %in% species_in_panel, is.finite(.data[[value_col]])) %>%
    mutate(plot_value = transform_metric_values(.data[[value_col]], transform_mode)) %>%
    filter(is.finite(.data$plot_value)) %>%
    add_numeric_positions(species_in_panel, geom_mode = geom_mode)

  local_scale <- identical(scale_mode, "species_local") && length(species_in_panel) > 1

  if (local_scale) {
    local_obj <- build_local_axis_df(
      plot_df,
      species_in_panel,
      transform_mode = transform_mode,
      ratio_metric = ratio_metric,
      geom_mode = geom_mode,
      negative_pad = use_negative_pad,
      axis_power = axis_power
    )
    plot_df <- plot_df %>%
      left_join(local_obj$win_df, by = c("species_chr")) %>%
      mutate(
        plot_value_clipped = pmin(pmax(.data$plot_value, .data$lower), .data$upper_core),
        plot_value_display = ifelse(
          .data$upper_core > .data$lower,
          LOCAL_BAND_BOTTOM + (.data$plot_value_clipped - .data$lower) / (.data$upper_core - .data$lower) * (LOCAL_BAND_TOP - LOCAL_BAND_BOTTOM),
          (LOCAL_BAND_BOTTOM + LOCAL_BAND_TOP) / 2
        )
      )
    sig_plot <- make_sig_numeric(stat_tbl, species_in_panel, LOCAL_SIG_BASE, LOCAL_SIG_STEP, local_scale = TRUE, geom_mode = geom_mode)
    n_df <- make_n_labels(plot_df, 0, 1, local_scale = TRUE)
    axis_df <- local_obj$axis_df
    plot_y <- "plot_value_display"
    local_lower <- max(0, LOCAL_BAND_BOTTOM - LOCAL_BOTTOM_MARGIN)
    local_upper <- if (nrow(sig_plot) > 0) {
      max(sig_plot$y + LOCAL_TOP_MARGIN, na.rm = TRUE)
    } else {
      LOCAL_BAND_TOP + LOCAL_TOP_MARGIN
    }
    if (!is.finite(local_upper)) local_upper <- LOCAL_BAND_TOP + LOCAL_TOP_MARGIN
    local_upper <- max(local_upper, LOCAL_BAND_TOP + LOCAL_TOP_MARGIN)
    local_upper <- min(local_upper, 0.94)
    y_limits <- c(local_lower, local_upper)
    y_expand_bottom <- 0
    y_expand_top <- 0
  } else {
    y_win <- get_display_window(plot_df$plot_value, transform_mode = transform_mode, ratio_metric = ratio_metric, negative_pad = use_negative_pad)
    sig_plot <- make_sig_numeric(stat_tbl, species_in_panel, y_win$base, y_win$step, local_scale = FALSE, geom_mode = geom_mode)
    y_upper <- if (nrow(sig_plot) > 0) {
      max(sig_plot$y + 1.2 * sig_plot$tip, na.rm = TRUE)
    } else {
      y_win$upper_core + 0.08 * y_win$range
    }
    if (!is.finite(y_upper)) y_upper <- y_win$upper
    y_upper <- max(y_upper, y_win$upper_core + 0.040 * y_win$range)
    n_df <- make_n_labels(plot_df, y_win$lower, y_win$range, local_scale = FALSE)
    axis_df <- NULL
    plot_df <- plot_df %>% mutate(plot_value_display = .data$plot_value)
    plot_y <- "plot_value_display"
    lower_pad <- ifelse(show_n && geom_mode == "box" && SHOW_N_IN_BOXPLOT, 0.07 * y_win$range, 0)
    y_limits <- c(y_win$lower - lower_pad, y_upper)
    y_expand_bottom <- 0
    y_expand_top <- 0.005
  }

  p <- ggplot(plot_df, aes(x = .data$xpos, y = .data[[plot_y]], fill = .data$fill_key, color = .data$fill_key))

  if (geom_mode == "violin_box") {
    p <- p +
      geom_violin(
        aes(group = interaction(.data$species, .data$gene_class)),
        width = cfg$geom_width,
        scale = "width",
        trim = TRUE,
        linewidth = 0.22,
        alpha = 0.90
      ) +
      geom_boxplot(
        aes(group = interaction(.data$species, .data$gene_class)),
        width = cfg$inner_width,
        coef = 0,
        fill = "white",
        color = "black",
        outlier.shape = NA,
        linewidth = 0.30
      )
  } else if (geom_mode == "box") {
    p <- p +
      geom_boxplot(
        aes(group = interaction(.data$species, .data$gene_class)),
        width = cfg$geom_width,
        coef = 1.5,
        outlier.shape = NA,
        linewidth = 0.28,
        alpha = 0.88
      )
  } else {
    stop("Unknown geom_mode: ", geom_mode, call. = FALSE)
  }

  p <- p +
    stat_summary(
      aes(group = interaction(.data$species, .data$gene_class)),
      fun = median,
      geom = "point",
      shape = 21,
      size = 0.70,
      stroke = 0.16,
      fill = "black",
      color = "black"
    )

  if (nrow(sig_plot) > 0) {
    p <- p +
      geom_segment(
        data = sig_plot,
        aes(x = .data$xmin, xend = .data$xmax, y = .data$y, yend = .data$y),
        inherit.aes = FALSE,
        linewidth = 0.22,
        color = "black"
      ) +
      geom_segment(
        data = sig_plot,
        aes(x = .data$xmin, xend = .data$xmin, y = .data$y, yend = .data$y - .data$tip),
        inherit.aes = FALSE,
        linewidth = 0.22,
        color = "black"
      ) +
      geom_segment(
        data = sig_plot,
        aes(x = .data$xmax, xend = .data$xmax, y = .data$y, yend = .data$y - .data$tip),
        inherit.aes = FALSE,
        linewidth = 0.22,
        color = "black"
      ) +
      geom_text(
        data = sig_plot,
        aes(x = (.data$xmin + .data$xmax) / 2, y = .data$y + 0.18 * .data$tip, label = .data$label),
        inherit.aes = FALSE,
        size = 2.05,
        color = "black",
        vjust = 0
      )
  }

  x_breaks <- unname(get_species_centers(species_in_panel, geom_mode))
  x_labels <- if (show_x_species_labels) parse(text = unname(species_axis_display[species_in_panel])) else rep("", length(species_in_panel))
  x_pad <- cfg$class_step + cfg$geom_width * 1.85
  x_min <- min(x_breaks, na.rm = TRUE) - x_pad
  x_max <- max(x_breaks, na.rm = TRUE) + x_pad

  if (!is.null(axis_df) && nrow(axis_df) > 0) {
    sep_positions <- if (length(x_breaks) >= 2) head(x_breaks, -1) + diff(x_breaks) / 2 else numeric(0)
    axis_position_map <- NULL
    axis_dir_map <- NULL
    if (length(species_in_panel) == 1) {
      axis_position_map <- setNames(x_min, species_in_panel[1])
      axis_dir_map <- setNames(-1, species_in_panel[1])
    } else if (length(species_in_panel) == 2) {
      axis_position_map <- setNames(c(x_min, x_max), species_in_panel[c(1, 2)])
      axis_dir_map <- setNames(c(-1, 1), species_in_panel[c(1, 2)])
    } else {
      axis_position_map <- setNames(c(x_min, sep_positions[1], x_max), species_in_panel[c(1, 2, 3)])
      axis_dir_map <- setNames(c(-1, -1, 1), species_in_panel[c(1, 2, 3)])
    }
    axis_df <- axis_df %>%
      mutate(
        axis_x = unname(axis_position_map[.data$species_chr]),
        axis_dir = unname(axis_dir_map[.data$species_chr]),
        tick_x = .data$axis_x + .data$axis_dir * cfg$geom_width * 0.26,
        label_x = .data$axis_x + .data$axis_dir * cfg$geom_width * 0.36,
        label_hjust = ifelse(.data$axis_dir < 0, 1, 0)
      )

    ## Full-height Y-axis spines: each species axis line runs across the whole
    ## panel (local_lower .. local_upper), so the left border, the first
    ## separator, and the right border read as real Y axes rather than short
    ## local segments. No separate separator segments are drawn, so exactly
    ## these three spines appear.
    axis_line_df <- axis_df %>%
      distinct(.data$species_chr, .data$axis_x) %>%
      mutate(ymin = local_lower, ymax = local_upper)

    p <- p +
      geom_segment(
        data = axis_line_df,
        aes(x = .data$axis_x, xend = .data$axis_x, y = .data$ymin, yend = .data$ymax),
        inherit.aes = FALSE,
        linewidth = 0.26,
        color = "black"
      ) +
      geom_segment(
        data = axis_df,
        aes(x = .data$axis_x, xend = .data$tick_x, y = .data$y_scaled, yend = .data$y_scaled),
        inherit.aes = FALSE,
        linewidth = 0.20,
        color = "black"
      ) +
      geom_text(
        data = axis_df,
        aes(x = .data$label_x, y = .data$y_scaled, label = .data$tick_label, hjust = .data$label_hjust),
        inherit.aes = FALSE,
        vjust = 0.5,
        size = 1.72,
        color = "black"
      )
  } else {
    x_breaks <- unname(get_species_centers(species_in_panel, geom_mode))
    x_labels <- if (show_x_species_labels) parse(text = unname(species_axis_display[species_in_panel])) else rep("", length(species_in_panel))
    x_pad <- cfg$class_step + cfg$geom_width * 1.85
    x_min <- min(x_breaks, na.rm = TRUE) - x_pad
    x_max <- max(x_breaks, na.rm = TRUE) + x_pad
  }

  if (show_n && SHOW_N_IN_BOXPLOT) {
    p <- p +
      geom_text(
        data = n_df,
        aes(x = .data$xpos, y = .data$y, label = .data$label),
        inherit.aes = FALSE,
        size = 1.70,
        angle = 0,
        hjust = 0.5,
        vjust = 1,
        color = "#333333"
      )
  }

  p <- p +
    scale_fill_manual(values = all_fill_values, guide = "none") +
    scale_color_manual(values = all_line_values, guide = "none") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(x_min, x_max),
      expand = expansion(mult = c(0.00, 0.00))
    ) +
    labs(x = NULL, y = ylab, title = title) +
    theme_compact(base_size = base_size)

  if (local_scale) {
    p <- p +
      scale_y_continuous(breaks = NULL, labels = NULL, limits = y_limits, expand = expansion(mult = c(0.00, 0.00))) +
      coord_cartesian(clip = "off")
  } else {
    p <- p +
      scale_y_continuous(
        breaks = function(x) nice_axis_breaks(x, axis_power = axis_power, n = 4),
        labels = function(x) fmt_axis_num(x, axis_power = axis_power),
        expand = expansion(mult = c(y_expand_bottom, y_expand_top))
      ) +
      coord_cartesian(ylim = y_limits, clip = "off")
  }

  p
}

class_legend_plot <- function(base_size = 7.4) {
  df <- data.frame(gene_class = factor(class_order, levels = class_order), x = 1, y = 1)
  ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$gene_class)) +
    geom_point(shape = 22, size = 3.0, color = "black", stroke = 0.20) +
    scale_fill_manual(
      values = class_legend_fill,
      breaks = class_order,
      name = "Gene class shade"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) +
    theme_void(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 0.2, color = "black"),
      legend.text = element_text(size = base_size - 0.3, color = "black"),
      legend.key.size = unit(0.32, "cm"),
      legend.spacing.x = unit(0.16, "cm"),
      plot.margin = margin(0, 0, 0, 0)
    )
}

combine_row_shared_y <- function(plot_list, y_title, title = NULL, tag_levels = NULL) {
  if (length(plot_list) == 0) stop("No plots supplied.", call. = FALSE)
  plot_list <- lapply(seq_along(plot_list), function(i) {
    if (i == 1) plot_list[[i]] + labs(y = y_title)
    else plot_list[[i]] + labs(y = NULL) + theme(axis.title.y = element_blank())
  })
  out <- wrap_plots(plot_list, nrow = 1)
  if (!is.null(title)) out <- out + plot_annotation(title = title)
  if (!is.null(tag_levels)) out <- out + plot_annotation(tag_levels = tag_levels)
  out
}

center_single_panel_row <- function(p, left = 1, center = 1, right = 1) {
  wrap_plots(plot_spacer(), p, plot_spacer(), nrow = 1, widths = c(left, center, right))
}

make_species_snp_row <- function(snp_df, snp_stats, data_type = PLOT_DATA_TYPE,
                                 geom_mode = "box", transform_mode = "raw") {
  df_plot <- snp_df %>%
    filter(as.character(.data$threshold) == SNP_DENSITY_THRESHOLD, .data$data_type == data_type)
  stat_plot <- snp_stats %>%
    filter(.data$threshold == SNP_DENSITY_THRESHOLD, .data$data_type == data_type, .data$value_col == "cds_snp_density")

  plots <- lapply(species_order, function(sp_id) {
    plot_numeric_panel(
      df = df_plot %>% filter(as.character(.data$species) == sp_id),
      stat_tbl = stat_plot %>% filter(.data$species == sp_id),
      value_col = "cds_snp_density",
      title = get_species_label(sp_id),
      ylab = NULL,
      species_in_panel = sp_id,
      transform_mode = transform_mode,
      geom_mode = geom_mode,
      show_x_species_labels = TRUE,
      show_n = geom_mode == "box",
      base_size = 7.4,
      ratio_metric = FALSE,
      scale_mode = "global"
    )
  })

  ylab <- make_y_label("CDS SNP density", "cds_snp_density", transform_mode)
  combine_row_shared_y(plots, ylab)
}

make_snp_combined_panel <- function(snp_df, snp_stats, data_type = PLOT_DATA_TYPE,
                                    geom_mode = "box", transform_mode = "raw",
                                    show_n = FALSE, scale_mode = "species_local") {
  df_plot <- snp_df %>%
    filter(as.character(.data$threshold) == SNP_DENSITY_THRESHOLD, .data$data_type == data_type)
  stat_plot <- snp_stats %>%
    filter(.data$threshold == SNP_DENSITY_THRESHOLD, .data$data_type == data_type, .data$value_col == "cds_snp_density")

  ylab <- make_y_label("CDS SNP density", "cds_snp_density", transform_mode)

  plot_numeric_panel(
    df = df_plot,
    stat_tbl = stat_plot,
    value_col = "cds_snp_density",
    title = paste0("SNP density, GERP > ", sub("gt", "", SNP_DENSITY_THRESHOLD)),
    ylab = ylab,
    species_in_panel = species_order,
    transform_mode = transform_mode,
    geom_mode = geom_mode,
    show_x_species_labels = TRUE,
    show_n = show_n,
    base_size = 7.4,
    ratio_metric = FALSE,
    scale_mode = scale_mode
  )
}

make_threshold_row <- function(df, stat_tbl, value_col, metric_title, ylab,
                               geom_mode = "box", transform_mode = "raw", ratio_metric = FALSE,
                               scale_mode = "species_local") {
  plots <- lapply(THRESHOLDS, function(thr) {
    plot_numeric_panel(
      df = df %>% filter(as.character(.data$threshold) == thr),
      stat_tbl = stat_tbl %>% filter(.data$threshold == thr, .data$value_col == value_col),
      value_col = value_col,
      title = paste0(metric_title, ", GERP > ", sub("gt", "", thr)),
      ylab = NULL,
      species_in_panel = species_order,
      transform_mode = transform_mode,
      geom_mode = geom_mode,
      show_x_species_labels = TRUE,
      show_n = geom_mode == "box",
      base_size = 7.4,
      ratio_metric = ratio_metric,
      scale_mode = scale_mode
    )
  })

  combine_row_shared_y(plots, ylab)
}

make_single_species_threshold_row <- function(df, stat_tbl, value_col, metric_title, ylab,
                                              sp_id,
                                              geom_mode = "box", transform_mode = "raw", ratio_metric = FALSE) {
  plots <- lapply(THRESHOLDS, function(thr) {
    plot_numeric_panel(
      df = df %>% filter(as.character(.data$threshold) == thr, as.character(.data$species) == sp_id),
      stat_tbl = stat_tbl %>% filter(.data$threshold == thr, .data$value_col == value_col, .data$species == sp_id),
      value_col = value_col,
      title = paste0(metric_title, ", GERP > ", sub("gt", "", thr)),
      ylab = NULL,
      species_in_panel = sp_id,
      transform_mode = transform_mode,
      geom_mode = geom_mode,
      show_x_species_labels = TRUE,
      show_n = geom_mode == "box",
      base_size = 7.4,
      ratio_metric = ratio_metric,
      scale_mode = "global"
    )
  })
  combine_row_shared_y(plots, ylab)
}

make_single_species_snp_panel <- function(snp_df, snp_stats, sp_id, data_type = PLOT_DATA_TYPE,
                                          geom_mode = "box", transform_mode = "raw") {
  df_plot <- snp_df %>% filter(as.character(.data$threshold) == SNP_DENSITY_THRESHOLD,
                               .data$data_type == data_type,
                               as.character(.data$species) == sp_id)
  stat_plot <- snp_stats %>% filter(.data$threshold == SNP_DENSITY_THRESHOLD,
                                    .data$data_type == data_type,
                                    .data$value_col == "cds_snp_density",
                                    .data$species == sp_id)
  ylab <- make_y_label("CDS SNP density", "cds_snp_density", transform_mode)
  plot_numeric_panel(
    df = df_plot,
    stat_tbl = stat_plot,
    value_col = "cds_snp_density",
    title = paste0(get_species_label(sp_id), ": SNP density, GERP > ", sub("gt", "", SNP_DENSITY_THRESHOLD)),
    ylab = ylab,
    species_in_panel = sp_id,
    transform_mode = transform_mode,
    geom_mode = geom_mode,
    show_x_species_labels = TRUE,
    show_n = geom_mode == "box",
    base_size = 7.4,
    ratio_metric = FALSE,
    scale_mode = "global"
  )
}

make_single_species_figure <- function(sp_id, geom_mode, mode,
                                       snp_metrics_all, snp_stats_all,
                                       gerp_metrics_all, gerp_stats_all,
                                       data_type = PLOT_DATA_TYPE) {
  ylab_gerp <- make_y_label("Constrained CDS sites (%)", "cds_conserved_pct", mode)
  ylab_den <- make_y_label("Deleterious SNP density", "cds_deleterious_density", mode)
  ylab_ratio <- make_y_label("Deleterious SNP ratio", "cds_deleterious_ratio", mode)

  p_gerp <- make_single_species_threshold_row(
    df = gerp_metrics_all,
    stat_tbl = gerp_stats_all,
    value_col = "cds_conserved_pct",
    metric_title = "Constrained CDS sites",
    ylab = ylab_gerp,
    sp_id = sp_id,
    geom_mode = geom_mode,
    transform_mode = mode,
    ratio_metric = FALSE
  )

  p_snp <- make_single_species_snp_panel(
    snp_df = snp_metrics_all,
    snp_stats = snp_stats_all,
    sp_id = sp_id,
    data_type = data_type,
    geom_mode = geom_mode,
    transform_mode = mode
  )

  p_den <- make_single_species_threshold_row(
    df = snp_metrics_all %>% filter(.data$data_type == data_type),
    stat_tbl = snp_stats_all %>% filter(.data$data_type == data_type),
    value_col = "cds_deleterious_density",
    metric_title = "Deleterious density",
    ylab = ylab_den,
    sp_id = sp_id,
    geom_mode = geom_mode,
    transform_mode = mode,
    ratio_metric = FALSE
  )

  p_ratio <- make_single_species_threshold_row(
    df = snp_metrics_all %>% filter(.data$data_type == data_type),
    stat_tbl = snp_stats_all %>% filter(.data$data_type == data_type),
    value_col = "cds_deleterious_ratio",
    metric_title = "Deleterious ratio",
    ylab = ylab_ratio,
    sp_id = sp_id,
    geom_mode = geom_mode,
    transform_mode = mode,
    ratio_metric = mode == "raw"
  )

  legend_p <- class_legend_plot(base_size = 7.4)

  ((p_gerp / p_snp / p_den / p_ratio / legend_p) +
      plot_layout(heights = c(1, 1, 1, 1, 0.12)) +
      plot_annotation(
        title = get_species_label(sp_id),
        tag_levels = "A",
        theme = theme(
          plot.title = element_text(size = 10.5, hjust = 0.5, face = "bold"),
          plot.tag = element_text(size = 10, face = "bold")
        )
      ))
}


## =========================================================
## 8. Main workflow
## =========================================================

cat("[INFO] Using PAN file: ", PAN_FILE, "\n", sep = "")

map_list_snp <- list(
  ped = read_core_sg_map_ped(),
  Zmay = build_gene_class_map_nonped("Zmay"),
  Atau = build_gene_class_map_nonped("Atau")
)

map_list_gerp <- list(
  ped = build_gene_class_map_ped_genes_for_gerp(),
  Zmay = map_list_snp$Zmay,
  Atau = map_list_snp$Atau
)

## Build SNP/deleterious metric tables.
snp_all_by_dt <- list()
snp_stats_by_dt <- list()
snp_summary_by_dt <- list()

for (thr in THRESHOLDS) {
  df_thr <- bind_rows(lapply(species_order, function(sp_id) {
    build_snp_one_threshold_species(sp_id, thr, map_list_snp[[sp_id]])
  }))

  split_list <- split_by_snp_class(df_thr)

  for (dt in DATA_TYPES) {
    df_dt <- split_list[[dt]] %>% mutate(data_type = dt)
    key <- paste(thr, dt, sep = "_")

    snp_all_by_dt[[key]] <- df_dt
    snp_stats_by_dt[[key]] <- bind_rows(
      run_stats_metric(df_dt, "cds_snp_density", "cds_snp_density", dt),
      run_stats_metric(df_dt, "cds_deleterious_density", "cds_deleterious_density", dt),
      run_stats_metric(df_dt, "cds_deleterious_ratio", "cds_deleterious_ratio", dt)
    )
    snp_summary_by_dt[[key]] <- summarise_metric_table(df_dt, dt)

    write_tsv_base(df_dt, file.path(TABLE_DIR, paste0("SNP_metrics_", thr, "_", dt, ".tsv")))
    write_tsv_base(snp_stats_by_dt[[key]], file.path(STAT_DIR, paste0("SNP_stats_", thr, "_", dt, ".tsv")))
    write_tsv_base(snp_summary_by_dt[[key]], file.path(TABLE_DIR, paste0("SNP_summary_", thr, "_", dt, ".tsv")))

    cat("[OK] SNP ", thr, " ", dt, ": n_units=", nrow(df_dt), "\n", sep = "")
  }
}

snp_metrics_all <- bind_rows(snp_all_by_dt)
snp_stats_all <- bind_rows(snp_stats_by_dt)
snp_summary_all <- bind_rows(snp_summary_by_dt)

write_tsv_base(snp_metrics_all, file.path(TABLE_DIR, "SNP_metrics_all_species_thresholds_data_types.tsv"))
write_tsv_base(snp_stats_all, file.path(STAT_DIR, "SNP_stats_all_species_thresholds_data_types.tsv"))
write_tsv_base(snp_summary_all, file.path(TABLE_DIR, "SNP_summary_all_species_thresholds_data_types.tsv"))

## Build GERP constrained CDS percentage tables.
gerp_metrics_all <- bind_rows(lapply(THRESHOLDS, function(thr) {
  bind_rows(lapply(species_order, function(sp_id) {
    build_gerp_one_threshold_species(sp_id, thr, map_list_gerp[[sp_id]])
  }))
})) %>%
  mutate(data_type = "GERP_CDS")

gerp_stats_all <- run_stats_metric(gerp_metrics_all, "cds_conserved_pct", "cds_conserved_pct", "GERP_CDS")
gerp_summary_all <- summarise_metric_table(gerp_metrics_all, "GERP_CDS")

write_tsv_base(gerp_metrics_all, file.path(TABLE_DIR, "GERP_CDS_constrained_pct_all_species_thresholds.tsv"))
write_tsv_base(gerp_stats_all, file.path(STAT_DIR, "GERP_CDS_constrained_pct_stats_all_species_thresholds.tsv"))
write_tsv_base(gerp_summary_all, file.path(TABLE_DIR, "GERP_CDS_constrained_pct_summary_all_species_thresholds.tsv"))


## =========================================================
## 9. Save figures
## =========================================================

legend_p <- class_legend_plot(base_size = 7.4)

for (geom_mode in GEOM_MODES_TO_SAVE) {
  for (mode in c("raw", "log10")) {
    mode_tag <- ifelse(mode == "raw", "raw_q99", "log10_q01_q99")

    ## SNP density: output both split-species row and combined local-scale panel.
    for (dt in DATA_TYPES) {
      p_snp_split <- make_species_snp_row(
        snp_metrics_all,
        snp_stats_all,
        data_type = dt,
        geom_mode = geom_mode,
        transform_mode = mode
      )
      p_snp_split_out <- p_snp_split / legend_p + plot_layout(heights = c(1, 0.12))
      prefix <- paste0("Fig_SNP_density_split_species_", dt, "_", mode_tag, "_", geom_mode)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_snp_split_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_snp_split_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)

      p_snp_combined <- make_snp_combined_panel(
        snp_metrics_all,
        snp_stats_all,
        data_type = dt,
        geom_mode = geom_mode,
        transform_mode = mode,
        show_n = geom_mode == "box",
        scale_mode = "species_local"
      )
      p_snp_combined_out <- p_snp_combined / legend_p + plot_layout(heights = c(1, 0.12))
      prefix <- paste0("Fig_SNP_density_combined_species_", dt, "_", mode_tag, "_", geom_mode)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_snp_combined_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_snp_combined_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)

      ylab_den <- make_y_label("Deleterious SNP density", "cds_deleterious_density", mode)
      p_den <- make_threshold_row(
        df = snp_metrics_all %>% filter(.data$data_type == dt),
        stat_tbl = snp_stats_all %>% filter(.data$data_type == dt),
        value_col = "cds_deleterious_density",
        metric_title = "Deleterious density",
        ylab = ylab_den,
        geom_mode = geom_mode,
        transform_mode = mode,
        ratio_metric = FALSE,
        scale_mode = "species_local"
      )
      p_den_out <- p_den / legend_p + plot_layout(heights = c(1, 0.12))
      prefix <- paste0("Fig_deleterious_density_by_GERP_threshold_", dt, "_", mode_tag, "_", geom_mode)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_den_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_den_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)

      ylab_ratio <- make_y_label("Deleterious SNP ratio", "cds_deleterious_ratio", mode)
      p_ratio <- make_threshold_row(
        df = snp_metrics_all %>% filter(.data$data_type == dt),
        stat_tbl = snp_stats_all %>% filter(.data$data_type == dt),
        value_col = "cds_deleterious_ratio",
        metric_title = "Deleterious ratio",
        ylab = ylab_ratio,
        geom_mode = geom_mode,
        transform_mode = mode,
        ratio_metric = mode == "raw",
        scale_mode = "species_local"
      )
      p_ratio_out <- p_ratio / legend_p + plot_layout(heights = c(1, 0.12))
      prefix <- paste0("Fig_deleterious_ratio_by_GERP_threshold_", dt, "_", mode_tag, "_", geom_mode)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_ratio_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_ratio_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
    }

    ylab_gerp <- make_y_label("Constrained CDS sites (%)", "cds_conserved_pct", mode)
    p_gerp <- make_threshold_row(
      df = gerp_metrics_all,
      stat_tbl = gerp_stats_all,
      value_col = "cds_conserved_pct",
      metric_title = "Constrained CDS sites",
      ylab = ylab_gerp,
      geom_mode = geom_mode,
      transform_mode = mode,
      ratio_metric = FALSE,
      scale_mode = "species_local"
    )
    p_gerp_out <- p_gerp / legend_p + plot_layout(heights = c(1, 0.12))
    prefix <- paste0("Fig_GERP_CDS_constrained_pct_by_threshold_", mode_tag, "_", geom_mode)
    safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_gerp_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)
    safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_gerp_out, ROW_FIG_WIDTH, ROW_FIG_HEIGHT, FIG_DPI)

    ## Final A4 portrait figure.
    p_gerp_main <- make_threshold_row(
      df = gerp_metrics_all,
      stat_tbl = gerp_stats_all,
      value_col = "cds_conserved_pct",
      metric_title = "Constrained CDS sites",
      ylab = ylab_gerp,
      geom_mode = geom_mode,
      transform_mode = mode,
      ratio_metric = FALSE,
      scale_mode = "species_local"
    )

    p_snp_main <- make_snp_combined_panel(
      snp_metrics_all,
      snp_stats_all,
      data_type = PLOT_DATA_TYPE,
      geom_mode = geom_mode,
      transform_mode = mode,
      show_n = geom_mode == "box",
      scale_mode = "species_local"
    )
    p_snp_main_row <- center_single_panel_row(p_snp_main)

    ylab_den <- make_y_label("Deleterious SNP density", "cds_deleterious_density", mode)
    p_den_main <- make_threshold_row(
      df = snp_metrics_all %>% filter(.data$data_type == PLOT_DATA_TYPE),
      stat_tbl = snp_stats_all %>% filter(.data$data_type == PLOT_DATA_TYPE),
      value_col = "cds_deleterious_density",
      metric_title = "Deleterious density",
      ylab = ylab_den,
      geom_mode = geom_mode,
      transform_mode = mode,
      ratio_metric = FALSE,
      scale_mode = "species_local"
    )

    ylab_ratio <- make_y_label("Deleterious SNP ratio", "cds_deleterious_ratio", mode)
    p_ratio_main <- make_threshold_row(
      df = snp_metrics_all %>% filter(.data$data_type == PLOT_DATA_TYPE),
      stat_tbl = snp_stats_all %>% filter(.data$data_type == PLOT_DATA_TYPE),
      value_col = "cds_deleterious_ratio",
      metric_title = "Deleterious ratio",
      ylab = ylab_ratio,
      geom_mode = geom_mode,
      transform_mode = mode,
      ratio_metric = mode == "raw",
      scale_mode = "species_local"
    )

    final_a4 <- (
      p_gerp_main /
        p_snp_main_row /
        p_den_main /
        p_ratio_main /
        legend_p
    ) +
      plot_layout(heights = c(1, 1, 1, 1, 0.12)) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(
          plot.tag = element_text(size = 10, face = "bold"),
          plot.margin = margin(3.5, 8, 3.5, 8)
        )
      )

    prefix <- paste0("Fig_final_A4portrait_GERP_SNP_deleterious_density_ratio_", PLOT_DATA_TYPE, "_", mode_tag, "_", geom_mode)
    safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), final_a4, A4_WIDTH, A4_HEIGHT, FIG_DPI)
    safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), final_a4, A4_WIDTH, A4_HEIGHT, FIG_DPI)

    ## Output single-species full figures.
    for (sp_id in species_order) {
      p_one <- make_single_species_figure(
        sp_id = sp_id,
        geom_mode = geom_mode,
        mode = mode,
        snp_metrics_all = snp_metrics_all,
        snp_stats_all = snp_stats_all,
        gerp_metrics_all = gerp_metrics_all,
        gerp_stats_all = gerp_stats_all,
        data_type = PLOT_DATA_TYPE
      )
      prefix <- paste0("Fig_single_species_", sp_id, "_", PLOT_DATA_TYPE, "_", mode_tag, "_", geom_mode)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".pdf")), p_one, SINGLE_SPECIES_WIDTH, SINGLE_SPECIES_HEIGHT, FIG_DPI)
      safe_ggsave(file.path(FIG_DIR, paste0(prefix, ".png")), p_one, SINGLE_SPECIES_WIDTH, SINGLE_SPECIES_HEIGHT, FIG_DPI)
    }
  }
}

cat("[DONE] All tables, statistics, and figures written to: ", OUTDIR, "\n", sep = "")
