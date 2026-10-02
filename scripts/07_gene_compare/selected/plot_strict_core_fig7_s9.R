#!/usr/bin/env Rscript

## ============================================================================
## Strict-core cross-species analysis for moso bamboo, teosinte/maize,
## and Tausch's goatgrass
##
## Analysis unit
##   One strict 1:1:1:1 SG from core_5sp_1to1.SG.pan.
##   P. edulis is represented by the combined pedC + pedD homeolog pair;
##   Z. mays and A. tauschii are represented by one gene each.
##
## Metrics
##   1. Constrained CDS sites (%)             GERP > 2, 4, 6
##   2. CDS SNP density                       plotted once from gt4 input
##   3. Putatively deleterious SNP density    GERP > 2, 4, 6
##   4. Putatively deleterious SNP ratio      GERP > 2, 4, 6; SNP-positive only
##   5. SNP0 SG proportion                    separate paired binary analysis
##
## Statistics
##   Continuous metrics: Friedman test + paired Wilcoxon signed-rank tests;
##                       BH correction within metric/threshold/data type.
##   SNP0 proportions:   Cochran's Q + pairwise exact McNemar tests;
##                       BH correction within threshold.
##
## Display
##   The visualization style is transferred from the unified gene-class script:
##   compact numeric X positions, reference box/violin dimensions, species-local
##   Y axes in three-species combined panels, five outward ticks per species,
##   one overall significance bracket, q-based display clipping only, and
##   600-dpi PDF/PNG output. Statistics always use the full original values.
##
## Required packages
##   dplyr, tidyr, ggplot2, patchwork, grid
## ============================================================================

SCRIPT_VERSION <- "2026-07-14-strict-core-shared-y-whisker-v10"

get_running_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return("<interactive/source>")
  fp <- sub("^--file=", "", file_arg[1])
  normalizePath(fp, winslash = "/", mustWork = FALSE)
}

## Print useful context before loading packages. If no [SCRIPT] lines appear,
## the failure occurred in an R startup profile before this script was read.
options(warn = 1)
options(error = function() {
  cat("\n[ERROR] Unhandled R error. Traceback follows:\n", file = stderr())
  traceback(30)
  quit(save = "no", status = 1, runLast = FALSE)
})

message("[SCRIPT] Version: ", SCRIPT_VERSION)
message("[SCRIPT] File: ", get_running_script())
message("[SCRIPT] Working directory: ",
        normalizePath(getwd(), winslash = "/", mustWork = FALSE))
message("[SCRIPT] R version: ", R.version.string)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

## ============================================================================
## 1. User parameters
## ============================================================================

PAN_FILE_CANDIDATES <- c(
  "core_5sp_1to1.SG.pan",
  "../../../01.OG_gene/core_5sp_1to1.SG.pan"
)

THRESHOLDS <- c("gt2", "gt4", "gt6")
SNP_DENSITY_THRESHOLD <- "gt4"

## Scale of an existing cds_pct / cds_conserved_pct column:
##   percent  = already 0-100
##   fraction = 0-1, multiply by 100
##   auto     = infer using max <= 1; use only when the input scale is unknown
EXISTING_CDS_PCT_SCALE <- "percent"

OUTDIR    <- "strict_core_3species_SNP_GERP_shared_y_whisker_v10_results"
TABLE_DIR <- file.path(OUTDIR, "tables")
STAT_DIR  <- file.path(OUTDIR, "stats")
FIG_DIR   <- file.path(OUTDIR, "figures")
QC_DIR    <- file.path(TABLE_DIR, "input_QC")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR,    recursive = TRUE, showWarnings = FALSE)

GEOM_MODES_TO_SAVE <- c("box", "violin_box")
DENSITY_TRANSFORMS <- c("raw", "log10")

## Display-only clipping. Statistics and output tables always use full values.
## GERP uses a 1%-99% display window; strongly right-skewed density/ratio
## panels use a pooled upper q97.5 window so the main distributions are visible.
DISPLAY_QUANTILES_RAW <- list(
  cds_conserved_pct       = c(0.010, 0.990),
  cds_snp_density         = c(0.000, 0.975),
  cds_deleterious_density = c(0.000, 0.975),
  cds_deleterious_ratio   = c(0.000, 0.975)
)
DISPLAY_QUANTILES_LOG10 <- list(
  cds_conserved_pct       = c(0.010, 0.990),
  cds_snp_density         = c(0.010, 0.990),
  cds_deleterious_density = c(0.010, 0.990),
  cds_deleterious_ratio   = c(0.010, 0.990)
)

PSEUDO <- 1e-6
SHOW_NS <- FALSE
SHOW_N <- TRUE
SHOW_SIGNIFICANCE <- TRUE

## Only one overall significance bracket is drawn in each panel.
## Pairwise Wilcoxon/McNemar results are still retained in the statistics TSVs.
SIGNIFICANCE_DISPLAY <- "overall"
N_Y_BREAKS <- 5
SHARE_Y_ACROSS_THRESHOLDS <- FALSE

## Geometry dimensions. The box inside a violin is intentionally narrow.
VIOLIN_WIDTH <- 0.58
VIOLIN_BOX_WIDTH <- 0.065
BOX_ONLY_WIDTH <- 0.30
VIOLIN_LINEWIDTH <- 0.26
BOX_LINEWIDTH <- 0.29

FIG_DPI <- 600
ROW_WIDTH <- 7.20
ROW_HEIGHT <- 2.52
A4_WIDTH <- 8.27
A4_HEIGHT <- 11.69

species_order <- c("P.edulis", "Z.mays", "A.tauschii")

species_labels <- c(
  "P.edulis"   = "italic(P.~edulis)",
  "Z.mays"     = "italic(Z.~mays)",
  "A.tauschii" = "italic(A.~tauschii)"
)

## Identical to the Core shades in the unified Core/Syntenic/Nonsyntenic plot.
species_fill <- c(
  "P.edulis"   = "#2F74B8",
  "Z.mays"     = "#5F3B2E",
  "A.tauschii" = "#507D39"
)

species_line <- c(
  "P.edulis"   = "#174D83",
  "Z.mays"     = "#3F261E",
  "A.tauschii" = "#2F4E22"
)

## ============================================================================
## 2. Generic helpers
## ============================================================================

find_first_existing <- function(paths, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(
      "Cannot find ", label, ". Tried:\n  ",
      paste(paths, collapse = "\n  "),
      call. = FALSE
    )
  }
  hit[1]
}

PAN_FILE <- find_first_existing(PAN_FILE_CANDIDATES, "strict-core SG pan file")

read_tsv_base <- function(fp, header = TRUE) {
  if (!file.exists(fp)) stop("Missing input file: ", fp, call. = FALSE)
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
    fp,
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
        filename, plot = plot, width = width, height = height,
        units = "in", device = cairo_pdf, limitsize = FALSE
      )
      ok <- TRUE
    }, silent = TRUE)
    if (!ok) {
      ggsave(
        filename, plot = plot, width = width, height = height,
        units = "in", limitsize = FALSE
      )
    }
  } else {
    ggsave(
      filename, plot = plot, width = width, height = height,
      units = "in", dpi = dpi, limitsize = FALSE
    )
  }
}

strip_species_prefix <- function(x) {
  x <- as.character(x)
  x <- sub("^Atau_", "", x)
  x <- sub("^A\\.tauschii_", "", x)
  x <- sub("^Osat_", "", x)
  x <- sub("^osat_", "", x)
  x <- sub("^pedC_", "", x)
  x <- sub("^pedD_", "", x)
  x <- sub("^Zmay_", "", x)
  x <- sub("^Z\\.mays_", "", x)
  x
}

p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = prob, names = FALSE, na.rm = TRUE))
}

normalize_pct_0_100 <- function(x, fp, input_scale = EXISTING_CDS_PCT_SCALE) {
  input_scale <- match.arg(input_scale, c("percent", "fraction", "auto"))
  x <- suppressWarnings(as.numeric(x))
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) return(x)

  if (input_scale == "fraction") {
    x <- 100 * x
  } else if (input_scale == "auto" && max(finite_x, na.rm = TRUE) <= 1.000001) {
    x <- 100 * x
    message("[INFO] Converted 0-1 CDS fraction to 0-100 in: ", fp)
  }

  if (any(is.finite(x) & (x < -1e-8 | x > 100.000001))) {
    warning(
      "CDS constrained percentage falls outside 0-100 in ", fp,
      ". Check EXISTING_CDS_PCT_SCALE and upstream site counts.",
      call. = FALSE
    )
  }
  x
}

assert_unique_gene_ids <- function(df, fp, value_cols = NULL) {
  dup_ids <- unique(df$gene_id_norm[duplicated(df$gene_id_norm)])
  if (length(dup_ids) == 0) return(invisible(df))

  if (is.null(value_cols)) {
    value_cols <- setdiff(
      names(df), c("source_row", "gene_id_raw", "gene_id_norm")
    )
  }
  dup_df <- df %>% filter(.data$gene_id_norm %in% dup_ids)
  conflict <- dup_df %>%
    group_by(.data$gene_id_norm) %>%
    summarise(
      across(
        all_of(value_cols),
        ~ n_distinct(.x, na.rm = FALSE),
        .names = "nuniq_{.col}"
      ),
      .groups = "drop"
    ) %>%
    filter(if_any(starts_with("nuniq_"), ~ .x > 1))

  if (nrow(conflict) > 0) {
    stop(
      "Conflicting duplicated normalized gene IDs in ", fp,
      ". Aggregate transcript-level rows upstream or select one canonical ",
      "transcript. Examples: ",
      paste(head(conflict$gene_id_norm, 10), collapse = ", "),
      call. = FALSE
    )
  }

  warning(
    "Collapsed ", length(dup_ids),
    " duplicated gene IDs with identical values in ", fp, ".",
    call. = FALSE
  )
  invisible(df)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  sum(x)
}

sum2_complete <- function(x, y) {
  ifelse(is.finite(x) & is.finite(y), x + y, NA_real_)
}

or_false <- function(x, y) {
  dplyr::coalesce(as.logical(x), FALSE) |
    dplyr::coalesce(as.logical(y), FALSE)
}

make_status <- function(flag_df) {
  if (nrow(flag_df) == 0) return(character())
  apply(as.data.frame(flag_df), 1, function(z) {
    bad <- names(z)[as.logical(z)]
    if (length(bad) == 0) "valid" else paste(bad, collapse = ";")
  })
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

qc_env <- new.env(parent = emptyenv())
qc_env$snp_input_summary <- list()
qc_env$gerp_input_summary <- list()
qc_env$snp_join_summary <- list()
qc_env$gerp_join_summary <- list()

record_qc <- function(slot, key, value) {
  obj <- qc_env[[slot]]
  obj[[key]] <- value
  qc_env[[slot]] <- obj
  invisible(value)
}

flush_qc_summaries <- function() {
  slots <- c(
    snp_input_summary = "SNP_input_QC_summary.tsv",
    gerp_input_summary = "GERP_input_QC_summary.tsv",
    snp_join_summary = "SNP_strict_SG_join_completeness.tsv",
    gerp_join_summary = "GERP_strict_SG_join_completeness.tsv"
  )
  for (slot in names(slots)) {
    obj <- qc_env[[slot]]
    if (length(obj) > 0) {
      write_tsv_base(bind_rows(obj), file.path(QC_DIR, slots[[slot]]))
    }
  }
}

## ============================================================================
## 3. Input file builders and checks
## ============================================================================

snp_files <- function(threshold) {
  list(
    pedC = paste0("pedC_", threshold, "_core_gene_deleterious.tsv"),
    pedD = paste0("pedD_", threshold, "_core_gene_deleterious.tsv"),
    Zmay = paste0("Zmay_", threshold, "_core_gene_deleterious.tsv"),
    Atau = paste0("Atau_", threshold, "_core_gene_deleterious.tsv")
  )
}

gerp_files <- function(threshold) {
  list(
    pedC = paste0("pedC_gerp_", threshold, "_CDS_per_gene.tsv"),
    pedD = paste0("pedD_gerp_", threshold, "_CDS_per_gene.tsv"),
    Zmay = paste0("Zmay_gerp_", threshold, "_CDS_per_gene.tsv"),
    Atau = paste0("Atau_gerp_", threshold, "_CDS_per_gene.tsv")
  )
}

check_all_inputs <- function() {
  required <- c(
    PAN_FILE,
    unlist(lapply(THRESHOLDS, snp_files), use.names = FALSE),
    unlist(lapply(THRESHOLDS, gerp_files), use.names = FALSE)
  )
  required <- unique(required)
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(
      "Missing input files or broken symbolic links:\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
  message("[OK] All required input files exist.")
}

## ============================================================================
## 4. Strict SG map
## ============================================================================

read_sg_pan <- function(fp) {
  pan <- read_tsv_base(fp, header = FALSE)
  if (ncol(pan) < 9) {
    stop("SG pan file must contain at least 9 columns: ", fp, call. = FALSE)
  }

  pan <- pan[, 1:9]
  colnames(pan) <- c(
    "SG", "n_species", "n_genes", "flag",
    "Atau_gene", "Osat_gene", "pedC_gene", "pedD_gene", "Zmay_gene"
  )

  out <- pan %>%
    transmute(
      SG = as.character(.data$SG),
      Atau = strip_species_prefix(.data$Atau_gene),
      pedC = strip_species_prefix(.data$pedC_gene),
      pedD = strip_species_prefix(.data$pedD_gene),
      Zmay = strip_species_prefix(.data$Zmay_gene)
    ) %>%
    filter(
      !is.na(.data$SG), .data$SG != "",
      !is.na(.data$Atau), .data$Atau != "-", .data$Atau != "",
      !is.na(.data$pedC), .data$pedC != "-", .data$pedC != "",
      !is.na(.data$pedD), .data$pedD != "-", .data$pedD != "",
      !is.na(.data$Zmay), .data$Zmay != "-", .data$Zmay != ""
    )

  if (anyDuplicated(out$SG)) stop("Duplicated SG IDs in pan file.", call. = FALSE)
  for (cc in c("Atau", "pedC", "pedD", "Zmay")) {
    if (anyDuplicated(out[[cc]])) {
      stop("A gene maps to multiple strict SGs in column ", cc, ".", call. = FALSE)
    }
  }

  message("[OK] Loaded strict SGs: ", nrow(out))
  out
}

## ============================================================================
## 5. SNP and putatively deleterious SNP metrics
## ============================================================================

read_snp_gene_table <- function(fp, label, threshold) {
  df <- read_tsv_base(fp, header = TRUE)
  if (!"gene_id" %in% names(df)) names(df)[1] <- "gene_id"

  required <- c("gene_id", "cds_len", "cds_snp_total", "cds_deleterious_snp")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing columns in ", fp, ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  parsed <- tibble(
    source_row = seq_len(nrow(df)) + 1L,
    gene_id_raw = as.character(df$gene_id),
    gene_id_norm = strip_species_prefix(df$gene_id),
    cds_len_raw = suppressWarnings(as.numeric(df$cds_len)),
    cds_snp_total_raw = suppressWarnings(as.numeric(df$cds_snp_total)),
    cds_deleterious_snp_raw =
      suppressWarnings(as.numeric(df$cds_deleterious_snp))
  ) %>%
    mutate(
      missing_gene_id = is.na(.data$gene_id_norm) | .data$gene_id_norm == "",
      missing_cds_len = !is.finite(.data$cds_len_raw),
      zero_cds_len = is.finite(.data$cds_len_raw) & .data$cds_len_raw == 0,
      negative_cds_len = is.finite(.data$cds_len_raw) & .data$cds_len_raw < 0,
      missing_total_snp = !is.finite(.data$cds_snp_total_raw),
      negative_total_snp = is.finite(.data$cds_snp_total_raw) &
        .data$cds_snp_total_raw < 0,
      missing_deleterious_snp = !is.finite(.data$cds_deleterious_snp_raw),
      negative_deleterious_snp = is.finite(.data$cds_deleterious_snp_raw) &
        .data$cds_deleterious_snp_raw < 0,
      deleterious_gt_total = is.finite(.data$cds_deleterious_snp_raw) &
        is.finite(.data$cds_snp_total_raw) &
        .data$cds_deleterious_snp_raw > .data$cds_snp_total_raw
    )

  flag_cols <- c(
    "missing_gene_id", "missing_cds_len", "zero_cds_len", "negative_cds_len",
    "missing_total_snp", "negative_total_snp", "missing_deleterious_snp",
    "negative_deleterious_snp", "deleterious_gt_total"
  )
  parsed$status <- make_status(parsed[, flag_cols, drop = FALSE])

  invalid_rows <- parsed %>% filter(.data$status != "valid")
  if (nrow(invalid_rows) > 0) {
    invalid_fp <- file.path(
      QC_DIR,
      paste0("SNP_", sanitize_filename(label), "_", threshold,
             "_invalid_rows.tsv")
    )
    write_tsv_base(invalid_rows, invalid_fp)
    warning(
      "Input QC: ", nrow(invalid_rows),
      " invalid/non-analyzable rows in ", fp,
      ". They were written to ", invalid_fp,
      ". The workflow continues using valid fields and matched SGs.",
      call. = FALSE
    )
  }

  valid_gene <- parsed %>% filter(!.data$missing_gene_id)
  assert_unique_gene_ids(
    valid_gene, fp,
    value_cols = c(
      "cds_len_raw", "cds_snp_total_raw", "cds_deleterious_snp_raw"
    )
  )
  valid_gene <- valid_gene %>% distinct(.data$gene_id_norm, .keep_all = TRUE)

  summary <- tibble(
    threshold = threshold,
    source = label,
    file = fp,
    n_input_rows = nrow(parsed),
    n_unique_gene_ids = n_distinct(valid_gene$gene_id_norm),
    n_missing_gene_id = sum(parsed$missing_gene_id),
    n_missing_cds_len = sum(parsed$missing_cds_len),
    n_zero_cds_len = sum(parsed$zero_cds_len),
    n_negative_cds_len = sum(parsed$negative_cds_len),
    n_missing_total_snp = sum(parsed$missing_total_snp),
    n_negative_total_snp = sum(parsed$negative_total_snp),
    n_missing_deleterious_snp = sum(parsed$missing_deleterious_snp),
    n_negative_deleterious_snp = sum(parsed$negative_deleterious_snp),
    n_deleterious_gt_total = sum(parsed$deleterious_gt_total),
    n_fully_valid_rows = sum(parsed$status == "valid")
  )
  record_qc("snp_input_summary", paste(label, threshold, sep = "__"), summary)

  cleaned <- valid_gene %>%
    transmute(
      gene_id_norm = .data$gene_id_norm,
      present = TRUE,
      cds_len = ifelse(
        is.finite(.data$cds_len_raw) & .data$cds_len_raw > 0,
        .data$cds_len_raw, NA_real_
      ),
      cds_snp_total = ifelse(
        is.finite(.data$cds_snp_total_raw) & .data$cds_snp_total_raw >= 0,
        .data$cds_snp_total_raw, NA_real_
      ),
      cds_deleterious_snp = ifelse(
        is.finite(.data$cds_deleterious_snp_raw) &
          .data$cds_deleterious_snp_raw >= 0 &
          is.finite(.data$cds_snp_total_raw) &
          .data$cds_snp_total_raw >= 0 &
          .data$cds_deleterious_snp_raw <= .data$cds_snp_total_raw,
        .data$cds_deleterious_snp_raw, NA_real_
      ),
      qc_deleterious_gt_total = .data$deleterious_gt_total,
      qc_status = .data$status
    )

  names(cleaned)[-1] <- paste0(label, "_", names(cleaned)[-1])
  cleaned
}

build_snp_wide_one_threshold <- function(threshold, pan_map) {
  ff <- snp_files(threshold)
  pedC <- read_snp_gene_table(ff$pedC, "pedC", threshold)
  pedD <- read_snp_gene_table(ff$pedD, "pedD", threshold)
  Zmay <- read_snp_gene_table(ff$Zmay, "Zmay", threshold)
  Atau <- read_snp_gene_table(ff$Atau, "Atau", threshold)

  wide <- pan_map %>%
    left_join(pedC, by = c("pedC" = "gene_id_norm")) %>%
    left_join(pedD, by = c("pedD" = "gene_id_norm")) %>%
    left_join(Zmay, by = c("Zmay" = "gene_id_norm")) %>%
    left_join(Atau, by = c("Atau" = "gene_id_norm")) %>%
    mutate(
      P.edulis_source_present = coalesce(.data$pedC_present, FALSE) &
        coalesce(.data$pedD_present, FALSE),
      P.edulis_source_invalid_del_gt_total = or_false(
        .data$pedC_qc_deleterious_gt_total,
        .data$pedD_qc_deleterious_gt_total
      ),
      P.edulis_cds_len = sum2_complete(
        .data$pedC_cds_len, .data$pedD_cds_len
      ),
      P.edulis_cds_snp_total = sum2_complete(
        .data$pedC_cds_snp_total, .data$pedD_cds_snp_total
      ),
      P.edulis_cds_deleterious_snp = sum2_complete(
        .data$pedC_cds_deleterious_snp,
        .data$pedD_cds_deleterious_snp
      ),

      Z.mays_source_present = coalesce(.data$Zmay_present, FALSE),
      Z.mays_source_invalid_del_gt_total = coalesce(
        .data$Zmay_qc_deleterious_gt_total, FALSE
      ),
      Z.mays_cds_len = .data$Zmay_cds_len,
      Z.mays_cds_snp_total = .data$Zmay_cds_snp_total,
      Z.mays_cds_deleterious_snp = .data$Zmay_cds_deleterious_snp,

      A.tauschii_source_present = coalesce(.data$Atau_present, FALSE),
      A.tauschii_source_invalid_del_gt_total = coalesce(
        .data$Atau_qc_deleterious_gt_total, FALSE
      ),
      A.tauschii_cds_len = .data$Atau_cds_len,
      A.tauschii_cds_snp_total = .data$Atau_cds_snp_total,
      A.tauschii_cds_deleterious_snp = .data$Atau_cds_deleterious_snp,
      threshold = threshold
    )

  for (sp in species_order) {
    len_col <- paste0(sp, "_cds_len")
    snp_col <- paste0(sp, "_cds_snp_total")
    del_col <- paste0(sp, "_cds_deleterious_snp")
    source_invalid_col <- paste0(sp, "_source_invalid_del_gt_total")
    invalid_col <- paste0(sp, "_invalid_del_gt_total")
    valid_snp_col <- paste0(sp, "_valid_snp_record")
    valid_del_col <- paste0(sp, "_valid_deleterious_record")
    snp_den_col <- paste0(sp, "_cds_snp_density")
    del_den_col <- paste0(sp, "_cds_deleterious_density")
    ratio_col <- paste0(sp, "_cds_deleterious_ratio")

    wide[[invalid_col]] <- coalesce(wide[[source_invalid_col]], FALSE) |
      (is.finite(wide[[del_col]]) & is.finite(wide[[snp_col]]) &
         wide[[del_col]] > wide[[snp_col]])
    wide[[valid_snp_col]] <- is.finite(wide[[len_col]]) &
      wide[[len_col]] > 0 & is.finite(wide[[snp_col]]) &
      wide[[snp_col]] >= 0
    wide[[valid_del_col]] <- wide[[valid_snp_col]] &
      is.finite(wide[[del_col]]) & wide[[del_col]] >= 0 &
      wide[[del_col]] <= wide[[snp_col]] & !wide[[invalid_col]]

    wide[[snp_den_col]] <- ifelse(
      wide[[valid_snp_col]],
      wide[[snp_col]] / wide[[len_col]], NA_real_
    )
    wide[[del_den_col]] <- ifelse(
      wide[[valid_del_col]],
      wide[[del_col]] / wide[[len_col]], NA_real_
    )
    wide[[ratio_col]] <- ifelse(
      wide[[valid_del_col]] & wide[[snp_col]] > 0,
      wide[[del_col]] / wide[[snp_col]], NA_real_
    )
  }

  join_summary <- bind_rows(lapply(species_order, function(sp) {
    len <- wide[[paste0(sp, "_cds_len")]]
    snp <- wide[[paste0(sp, "_cds_snp_total")]]
    del <- wide[[paste0(sp, "_cds_deleterious_snp")]]
    present <- wide[[paste0(sp, "_source_present")]]
    source_bad <- wide[[paste0(sp, "_source_invalid_del_gt_total")]]
    valid_snp <- is.finite(len) & len > 0 & is.finite(snp) & snp >= 0
    valid_del <- valid_snp & is.finite(del) & del >= 0 &
      del <= snp & !source_bad
    tibble(
      threshold = threshold,
      species = sp,
      n_total_strict_SG = nrow(wide),
      n_source_present = sum(present, na.rm = TRUE),
      n_valid_SNP_density = sum(valid_snp),
      n_valid_deleterious_density = sum(valid_del),
      n_valid_ratio = sum(valid_del & snp > 0),
      n_valid_SNP0 = sum(valid_snp & snp == 0),
      n_source_deleterious_gt_total = sum(source_bad, na.rm = TRUE)
    )
  }))
  record_qc("snp_join_summary", threshold, join_summary)

  n_invalid <- sum(
    as.matrix(wide[, paste0(species_order, "_invalid_del_gt_total"),
                   drop = FALSE]),
    na.rm = TRUE
  )
  if (n_invalid > 0) {
    warning(
      "At ", threshold, ", detected ", n_invalid,
      " species-SG records with deleterious SNP count > total SNP count. ",
      "Their deleterious density and ratio are NA; inspect input_QC tables.",
      call. = FALSE
    )
  }

  wide
}

prepare_paired_metric_wide <- function(wide, metric, data_type = "ALL") {
  use <- wide
  snp_cols <- paste0(species_order, "_cds_snp_total")
  if (data_type == "SNPpos") {
    use <- use %>%
      filter(if_all(all_of(snp_cols), ~ is.finite(.x) & .x > 0))
  } else if (!data_type %in% c("ALL", "GERP_CDS")) {
    stop("Unknown data_type: ", data_type, call. = FALSE)
  }

  value_cols <- paste0(species_order, "_", metric)
  use %>%
    select(.data$SG, all_of(value_cols)) %>%
    rename_with(
      ~ sub(paste0("_", metric, "$"), "", .x),
      all_of(value_cols)
    ) %>%
    filter(if_all(all_of(species_order), ~ is.finite(.x)))
}

metric_long_from_wide <- function(wide, metric, data_type = "ALL") {
  paired <- prepare_paired_metric_wide(wide, metric, data_type)
  threshold_value <- unique(wide$threshold)
  threshold_value <- threshold_value[1]

  if (nrow(paired) == 0) {
    warning(
      "No complete matched SGs for metric=", metric,
      ", data_type=", data_type,
      ", threshold=", threshold_value, ".",
      call. = FALSE
    )
    return(tibble(
      SG = character(), threshold = character(),
      species = factor(character(), levels = species_order),
      value = numeric(), metric = character(), data_type = character()
    ))
  }

  paired %>%
    mutate(threshold = threshold_value) %>%
    pivot_longer(
      cols = all_of(species_order),
      names_to = "species", values_to = "value"
    ) %>%
    mutate(
      species = factor(.data$species, levels = species_order),
      metric = metric,
      data_type = data_type,
      value = as.numeric(.data$value)
    )
}

summarise_long_metric <- function(long_df) {
  long_df %>%
    group_by(.data$threshold, .data$data_type, .data$metric, .data$species) %>%
    summarise(
      n = n(),
      median = median(.data$value),
      q1 = safe_quantile(.data$value, 0.25),
      q3 = safe_quantile(.data$value, 0.75),
      mean = mean(.data$value),
      sd = sd(.data$value),
      min = min(.data$value),
      max = max(.data$value),
      .groups = "drop"
    )
}

ratio_audit <- function(wide) {
  bind_rows(lapply(species_order, function(sp) {
    snp <- wide[[paste0(sp, "_cds_snp_total")]]
    del <- wide[[paste0(sp, "_cds_deleterious_snp")]]
    ratio <- wide[[paste0(sp, "_cds_deleterious_ratio")]]
    invalid <- wide[[paste0(sp, "_invalid_del_gt_total")]]
    total_snp_sum <- safe_sum(snp)
    total_del_sum <- safe_sum(del)
    tibble(
      threshold = unique(wide$threshold)[1],
      species = sp,
      n_total_SG = nrow(wide),
      n_missing_or_invalid_total_SNP = sum(!is.finite(snp)),
      n_SNP0_SG = sum(snp == 0, na.rm = TRUE),
      n_SNP_positive_SG = sum(snp > 0, na.rm = TRUE),
      n_ratio_defined = sum(is.finite(ratio)),
      n_invalid_del_gt_total = sum(invalid, na.rm = TRUE),
      total_CDS_SNPs = total_snp_sum,
      total_deleterious_SNPs = total_del_sum,
      median_gene_ratio = safe_median(ratio),
      q1_gene_ratio = safe_quantile(ratio, 0.25),
      q3_gene_ratio = safe_quantile(ratio, 0.75),
      pooled_ratio = ifelse(
        is.finite(total_snp_sum) & total_snp_sum > 0 &
          is.finite(total_del_sum),
        total_del_sum / total_snp_sum, NA_real_
      )
    )
  }))
}

audit_snp_inputs_across_thresholds <- function(wide_list) {
  keys <- c(
    paste0(species_order, "_cds_len"),
    paste0(species_order, "_cds_snp_total")
  )
  ref <- wide_list[[SNP_DENSITY_THRESHOLD]] %>%
    select(.data$SG, all_of(keys))

  bind_rows(lapply(setdiff(THRESHOLDS, SNP_DENSITY_THRESHOLD), function(thr) {
    cur <- wide_list[[thr]] %>% select(.data$SG, all_of(keys))
    joined <- full_join(ref, cur, by = "SG", suffix = c("_ref", "_cur"))
    bind_rows(lapply(keys, function(k) {
      x <- joined[[paste0(k, "_ref")]]
      y <- joined[[paste0(k, "_cur")]]
      both <- is.finite(x) & is.finite(y)
      missing_pattern_diff <- xor(is.finite(x), is.finite(y))
      value_diff <- both & x != y
      tibble(
        reference_threshold = SNP_DENSITY_THRESHOLD,
        compared_threshold = thr,
        field = k,
        n_rows = nrow(joined),
        n_both_finite = sum(both),
        n_missing_reference = sum(!is.finite(x)),
        n_missing_compared = sum(!is.finite(y)),
        n_missing_pattern_different = sum(missing_pattern_diff),
        n_value_different = sum(value_diff),
        max_abs_difference = ifelse(
          any(value_diff),
          max(abs(x[value_diff] - y[value_diff]), na.rm = TRUE), 0
        )
      )
    }))
  }))
}

## ============================================================================
## 6. GERP constrained CDS percentage
## ============================================================================

read_gerp_gene_table <- function(fp, label, threshold) {
  df <- read_tsv_base(fp, header = TRUE)
  if (!"gene_id" %in% names(df)) names(df)[1] <- "gene_id"

  cds_len_raw <- if ("cds_len" %in% names(df)) {
    suppressWarnings(as.numeric(df$cds_len))
  } else {
    rep(NA_real_, nrow(df))
  }

  cds_sites_raw <- if ("cds_conserved_sites" %in% names(df)) {
    suppressWarnings(as.numeric(df$cds_conserved_sites))
  } else if ("cds_sites" %in% names(df)) {
    suppressWarnings(as.numeric(df$cds_sites))
  } else {
    rep(NA_real_, nrow(df))
  }

  cds_pct_existing <- if ("cds_conserved_pct" %in% names(df)) {
    normalize_pct_0_100(df$cds_conserved_pct, fp)
  } else if ("cds_pct" %in% names(df)) {
    normalize_pct_0_100(df$cds_pct, fp)
  } else {
    rep(NA_real_, nrow(df))
  }

  parsed <- tibble(
    source_row = seq_len(nrow(df)) + 1L,
    gene_id_raw = as.character(df$gene_id),
    gene_id_norm = strip_species_prefix(df$gene_id),
    cds_len_raw = cds_len_raw,
    cds_sites_raw = cds_sites_raw,
    cds_pct_existing = cds_pct_existing
  ) %>%
    mutate(
      missing_gene_id = is.na(.data$gene_id_norm) | .data$gene_id_norm == "",
      missing_cds_len = !is.finite(.data$cds_len_raw),
      zero_cds_len = is.finite(.data$cds_len_raw) & .data$cds_len_raw == 0,
      negative_cds_len = is.finite(.data$cds_len_raw) & .data$cds_len_raw < 0,
      negative_cds_sites = is.finite(.data$cds_sites_raw) &
        .data$cds_sites_raw < 0,
      cds_sites_gt_len = is.finite(.data$cds_sites_raw) &
        is.finite(.data$cds_len_raw) & .data$cds_len_raw >= 0 &
        .data$cds_sites_raw > .data$cds_len_raw,
      existing_pct_outside_0_100 = is.finite(.data$cds_pct_existing) &
        (.data$cds_pct_existing < 0 | .data$cds_pct_existing > 100),
      count_pct_valid = is.finite(.data$cds_sites_raw) &
        .data$cds_sites_raw >= 0 & is.finite(.data$cds_len_raw) &
        .data$cds_len_raw > 0 & .data$cds_sites_raw <= .data$cds_len_raw,
      existing_pct_valid = is.finite(.data$cds_pct_existing) &
        .data$cds_pct_existing >= 0 & .data$cds_pct_existing <= 100,
      cds_pct_clean = case_when(
        .data$count_pct_valid ~ 100 * .data$cds_sites_raw / .data$cds_len_raw,
        .data$existing_pct_valid ~ .data$cds_pct_existing,
        TRUE ~ NA_real_
      ),
      pct_source = case_when(
        .data$count_pct_valid ~ "site_count_over_CDS_length",
        .data$existing_pct_valid ~ "existing_percentage",
        TRUE ~ "unavailable"
      )
    )

  flag_cols <- c(
    "missing_gene_id", "negative_cds_len", "negative_cds_sites",
    "cds_sites_gt_len", "existing_pct_outside_0_100"
  )
  parsed$status <- make_status(parsed[, flag_cols, drop = FALSE])
  parsed$status[!is.finite(parsed$cds_pct_clean) &
                  parsed$status == "valid"] <- "no_valid_CDS_percentage"

  invalid_rows <- parsed %>%
    filter(.data$status != "valid" | !is.finite(.data$cds_pct_clean))
  if (nrow(invalid_rows) > 0) {
    invalid_fp <- file.path(
      QC_DIR,
      paste0("GERP_", sanitize_filename(label), "_", threshold,
             "_invalid_rows.tsv")
    )
    write_tsv_base(invalid_rows, invalid_fp)
    warning(
      "GERP input QC: ", nrow(invalid_rows),
      " invalid/non-analyzable rows in ", fp,
      ". They were written to ", invalid_fp, ".",
      call. = FALSE
    )
  }

  valid_gene <- parsed %>% filter(!.data$missing_gene_id)
  assert_unique_gene_ids(
    valid_gene, fp,
    value_cols = c("cds_len_raw", "cds_sites_raw", "cds_pct_existing")
  )
  valid_gene <- valid_gene %>% distinct(.data$gene_id_norm, .keep_all = TRUE)

  summary <- tibble(
    threshold = threshold,
    source = label,
    file = fp,
    n_input_rows = nrow(parsed),
    n_unique_gene_ids = n_distinct(valid_gene$gene_id_norm),
    n_missing_gene_id = sum(parsed$missing_gene_id),
    n_missing_cds_len = sum(parsed$missing_cds_len),
    n_zero_cds_len = sum(parsed$zero_cds_len),
    n_negative_cds_len = sum(parsed$negative_cds_len),
    n_negative_cds_sites = sum(parsed$negative_cds_sites),
    n_cds_sites_gt_len = sum(parsed$cds_sites_gt_len),
    n_existing_pct_outside_0_100 =
      sum(parsed$existing_pct_outside_0_100),
    n_valid_pct = sum(is.finite(parsed$cds_pct_clean)),
    n_pct_from_site_counts =
      sum(parsed$pct_source == "site_count_over_CDS_length"),
    n_pct_from_existing_column =
      sum(parsed$pct_source == "existing_percentage")
  )
  record_qc("gerp_input_summary", paste(label, threshold, sep = "__"), summary)

  cleaned <- valid_gene %>%
    transmute(
      gene_id_norm = .data$gene_id_norm,
      present = TRUE,
      cds_len = ifelse(
        is.finite(.data$cds_len_raw) & .data$cds_len_raw > 0,
        .data$cds_len_raw, NA_real_
      ),
      cds_sites = ifelse(
        .data$count_pct_valid, .data$cds_sites_raw, NA_real_
      ),
      cds_pct = .data$cds_pct_clean,
      pct_source = .data$pct_source,
      qc_status = .data$status
    )
  names(cleaned)[-1] <- paste0(label, "_", names(cleaned)[-1])
  cleaned
}

combine_ped_gerp_pct <- function(wide) {
  count_ok <- is.finite(wide$pedC_cds_sites) &
    is.finite(wide$pedD_cds_sites) & is.finite(wide$pedC_cds_len) &
    is.finite(wide$pedD_cds_len) &
    (wide$pedC_cds_len + wide$pedD_cds_len > 0)
  weighted_pct_ok <- is.finite(wide$pedC_cds_pct) &
    is.finite(wide$pedD_cds_pct) & is.finite(wide$pedC_cds_len) &
    is.finite(wide$pedD_cds_len) &
    (wide$pedC_cds_len + wide$pedD_cds_len > 0)
  mean_pct_ok <- is.finite(wide$pedC_cds_pct) &
    is.finite(wide$pedD_cds_pct)

  count_pct <- 100 *
    (wide$pedC_cds_sites + wide$pedD_cds_sites) /
    (wide$pedC_cds_len + wide$pedD_cds_len)
  weighted_pct <-
    (wide$pedC_cds_pct * wide$pedC_cds_len +
       wide$pedD_cds_pct * wide$pedD_cds_len) /
    (wide$pedC_cds_len + wide$pedD_cds_len)
  mean_pct <- rowMeans(
    wide[, c("pedC_cds_pct", "pedD_cds_pct")], na.rm = FALSE
  )

  tibble(
    value = ifelse(
      count_ok, count_pct,
      ifelse(
        weighted_pct_ok, weighted_pct,
        ifelse(mean_pct_ok, mean_pct, NA_real_)
      )
    ),
    method = ifelse(
      count_ok, "summed_sites_over_summed_length",
      ifelse(
        weighted_pct_ok, "CDS_length_weighted_pct",
        ifelse(mean_pct_ok, "unweighted_mean_pedC_pedD_pct", "unavailable")
      )
    )
  )
}

build_gerp_wide_one_threshold <- function(threshold, pan_map) {
  ff <- gerp_files(threshold)
  pedC <- read_gerp_gene_table(ff$pedC, "pedC", threshold)
  pedD <- read_gerp_gene_table(ff$pedD, "pedD", threshold)
  Zmay <- read_gerp_gene_table(ff$Zmay, "Zmay", threshold)
  Atau <- read_gerp_gene_table(ff$Atau, "Atau", threshold)

  wide <- pan_map %>%
    left_join(pedC, by = c("pedC" = "gene_id_norm")) %>%
    left_join(pedD, by = c("pedD" = "gene_id_norm")) %>%
    left_join(Zmay, by = c("Zmay" = "gene_id_norm")) %>%
    left_join(Atau, by = c("Atau" = "gene_id_norm"))

  ped <- combine_ped_gerp_pct(wide)
  wide <- wide %>%
    mutate(
      threshold = threshold,
      P.edulis_source_present = coalesce(.data$pedC_present, FALSE) &
        coalesce(.data$pedD_present, FALSE),
      P.edulis_cds_conserved_pct = ped$value,
      P.edulis_aggregation_method = ped$method,
      Z.mays_source_present = coalesce(.data$Zmay_present, FALSE),
      Z.mays_cds_conserved_pct = .data$Zmay_cds_pct,
      A.tauschii_source_present = coalesce(.data$Atau_present, FALSE),
      A.tauschii_cds_conserved_pct = .data$Atau_cds_pct
    )

  join_summary <- bind_rows(lapply(species_order, function(sp) {
    pct <- wide[[paste0(sp, "_cds_conserved_pct")]]
    present <- wide[[paste0(sp, "_source_present")]]
    tibble(
      threshold = threshold,
      species = sp,
      n_total_strict_SG = nrow(wide),
      n_source_present = sum(present, na.rm = TRUE),
      n_valid_constrained_pct = sum(
        is.finite(pct) & pct >= 0 & pct <= 100
      )
    )
  }))
  record_qc("gerp_join_summary", threshold, join_summary)

  n_fallback <- sum(
    wide$P.edulis_aggregation_method ==
      "unweighted_mean_pedC_pedD_pct",
    na.rm = TRUE
  )
  if (n_fallback > 0) {
    warning(
      "At ", threshold, ", ", n_fallback,
      " P. edulis SGs used an unweighted mean because valid CDS lengths ",
      "were unavailable. Inspect P.edulis_aggregation_method.",
      call. = FALSE
    )
  }

  wide
}

gerp_long_from_wide <- function(wide) {
  metric_long_from_wide(wide, "cds_conserved_pct", "GERP_CDS")
}

## ============================================================================
## 7. Paired statistics
## ============================================================================

paired_continuous_stats <- function(wide, metric, threshold, data_type) {
  stat_df <- prepare_paired_metric_wide(wide, metric, data_type)

  if (nrow(stat_df) < 3) {
    return(tibble(
      threshold = threshold, data_type = data_type, metric = metric,
      test = NA_character_, Friedman_chisq = NA_real_,
      Friedman_df = NA_real_, Friedman_p = NA_real_,
      Comparison = NA_character_, group1 = NA_character_,
      group2 = NA_character_, n_pairs = nrow(stat_df),
      median_group1 = NA_real_, median_group2 = NA_real_,
      median_paired_difference = NA_real_, P.raw = NA_real_,
      P.adj = NA_real_, label = NA_character_
    ))
  }

  mat <- as.matrix(stat_df[, species_order, drop = FALSE])
  storage.mode(mat) <- "numeric"
  fried <- tryCatch(friedman.test(mat), error = function(e) NULL)

  out <- bind_rows(lapply(combn(species_order, 2, simplify = FALSE), function(cc) {
    g1 <- cc[1]
    g2 <- cc[2]
    x <- stat_df[[g1]]
    y <- stat_df[[g2]]
    wt <- tryCatch(
      suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)),
      error = function(e) NULL
    )
    tibble(
      threshold = threshold,
      data_type = data_type,
      metric = metric,
      test = "Friedman + paired Wilcoxon signed-rank; BH correction",
      Friedman_chisq = ifelse(
        is.null(fried), NA_real_, unname(fried$statistic)
      ),
      Friedman_df = ifelse(
        is.null(fried), NA_real_, unname(fried$parameter)
      ),
      Friedman_p = ifelse(
        is.null(fried), NA_real_, fried$p.value
      ),
      Comparison = paste0(g1, " - ", g2),
      group1 = g1,
      group2 = g2,
      n_pairs = length(x),
      median_group1 = median(x),
      median_group2 = median(y),
      median_paired_difference = median(x - y),
      P.raw = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }))

  out %>%
    mutate(
      P.adj = p.adjust(.data$P.raw, method = "BH"),
      label = vapply(.data$P.adj, p_to_star, character(1))
    )
}

cochran_q_test_base <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[complete.cases(mat), , drop = FALSE]
  n <- nrow(mat)
  k <- ncol(mat)
  if (n < 2 || k < 2) {
    return(list(statistic = NA_real_, df = k - 1,
                p.value = NA_real_, n = n))
  }
  Cj <- colSums(mat)
  Ri <- rowSums(mat)
  total <- sum(Cj)
  denominator <- k * total - sum(Ri^2)
  q <- if (denominator > 0) {
    (k - 1) * (k * sum(Cj^2) - total^2) / denominator
  } else {
    NA_real_
  }
  list(
    statistic = q,
    df = k - 1,
    p.value = ifelse(
      is.finite(q), pchisq(q, df = k - 1, lower.tail = FALSE), NA_real_
    ),
    n = n
  )
}

snp0_paired_stats <- function(wide, source_threshold) {
  valid_cols <- paste0(species_order, "_valid_snp_record")
  count_cols <- paste0(species_order, "_cds_snp_total")

  dat <- wide %>%
    filter(if_all(all_of(valid_cols), ~ .x %in% TRUE)) %>%
    select(.data$SG, all_of(count_cols)) %>%
    rename_with(~ sub("_cds_snp_total$", "", .x), all_of(count_cols))

  binary <- as.data.frame(lapply(
    dat[, species_order, drop = FALSE], function(x) as.integer(x == 0)
  ))
  q <- cochran_q_test_base(binary)

  pair_tbl <- bind_rows(lapply(
    combn(species_order, 2, simplify = FALSE), function(cc) {
      g1 <- cc[1]
      g2 <- cc[2]
      x <- binary[[g1]]
      y <- binary[[g2]]
      b <- sum(x == 1 & y == 0)
      c <- sum(x == 0 & y == 1)
      discordant <- b + c
      exact_p <- if (discordant > 0) {
        binom.test(
          b, discordant, p = 0.5, alternative = "two.sided"
        )$p.value
      } else {
        1
      }
      tibble(
        source_threshold = source_threshold,
        metric = "SNP0_SG_proportion",
        test = "Cochran Q + pairwise exact McNemar; BH correction",
        Cochran_Q = q$statistic,
        Cochran_df = q$df,
        Cochran_p = q$p.value,
        n_complete_SG = q$n,
        Comparison = paste0(g1, " - ", g2),
        group1 = g1,
        group2 = g2,
        n_group1_SNP0 = sum(x == 1),
        n_group2_SNP0 = sum(y == 1),
        prop_group1_SNP0 = mean(x),
        prop_group2_SNP0 = mean(y),
        discordant_group1_only = b,
        discordant_group2_only = c,
        discordant_odds_ratio = ifelse(
          c > 0, b / c, ifelse(b > 0, Inf, NA_real_)
        ),
        P.raw = exact_p
      )
    }
  )) %>%
    mutate(
      P.adj = p.adjust(.data$P.raw, method = "BH"),
      label = vapply(.data$P.adj, p_to_star, character(1))
    )

  summary <- tibble(
    source_threshold = source_threshold,
    species = factor(species_order, levels = species_order),
    n_complete_SG = nrow(binary),
    n_SNP0_SG = vapply(
      binary[, species_order, drop = FALSE], sum, numeric(1)
    ),
    n_SNP_positive_SG = nrow(binary) - vapply(
      binary[, species_order, drop = FALSE], sum, numeric(1)
    ),
    prop_SNP0 = vapply(
      binary[, species_order, drop = FALSE], mean, numeric(1)
    )
  )

  list(summary = summary, stats = pair_tbl)
}

## ============================================================================
## 8. Plot helpers: shared Y based on Tukey whiskers + fixed pairwise order
## ============================================================================

## This plotting layer follows the updated display rules:
##   1. one shared Y axis within each multi-threshold row;
##   2. boxplot Y-range is determined from the Tukey whiskers of the three
##      species combined by taking min(lower whisker) and max(upper whisker);
##   3. outliers beyond whiskers are hidden;
##   4. violin tails may be display-clipped (q-based), but the overlaid box and
##      whiskers always remain fully visible;
##   5. pairwise significance brackets are shown in a fixed bottom-to-top order:
##        ped vs Z.mays, ped vs A.tauschii, Z.mays vs A.tauschii.

CONTINUOUS_DATA_TYPES <- c("ALL", "SNPpos")

REF_BOX_WIDTH <- 0.34
REF_VIOLIN_WIDTH <- 0.58
REF_VIOLIN_BOX_WIDTH <- 0.13
REF_SNP0_BAR_WIDTH <- 0.34
REF_AXIS_TICK_LINEWIDTH <- 0.24
REF_AXIS_TICK_LENGTH_PT <- 1.35
REF_N_Y_BREAKS <- 5

ROW_WIDTH <- 7.60
ROW_HEIGHT <- 2.45
A4_WIDTH <- 8.27
A4_HEIGHT <- 11.69
FIG_DPI <- 600

SHOW_N <- TRUE
SHOW_SIGNIFICANCE <- TRUE
SHOW_NS <- FALSE
SIGNIFICANCE_DISPLAY <- "pairwise"

REF_DISPLAY_QUANTILES_RAW <- list(
  cds_conserved_pct = c(0.000, 0.990),
  cds_snp_density = c(0.000, 0.975),
  cds_deleterious_density = c(0.000, 0.975),
  cds_deleterious_ratio = c(0.000, 0.975)
)
REF_DISPLAY_QUANTILES_LOG10 <- list(
  cds_conserved_pct = c(0.010, 0.990),
  cds_snp_density = c(0.010, 0.990),
  cds_deleterious_density = c(0.010, 0.990),
  cds_deleterious_ratio = c(0.010, 0.990)
)

species_axis_display <- c(
  "P.edulis" = "italic(P.~edulis)",
  "Z.mays" = "italic(Z.~mays)",
  "A.tauschii" = "italic(A.~tauschii)"
)

pairwise_display_order_ref <- tribble(
  ~group1, ~group2, ~ord,
  "P.edulis", "Z.mays", 1,
  "P.edulis", "A.tauschii", 2,
  "Z.mays", "A.tauschii", 3
)

safe_quantile_ref <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

theme_core_reference <- function(base_size = 7.4) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 0.4, face = "plain", hjust = 0.5,
        margin = margin(b = 2.5)
      ),
      axis.text.x = element_text(
        size = base_size - 0.35, color = "black", margin = margin(t = 0.9)
      ),
      axis.text.y = element_text(size = base_size - 0.5, color = "black"),
      axis.title.y = element_text(
        size = base_size + 0.25, color = "black", margin = margin(r = 4.0)
      ),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.28),
      axis.ticks = element_line(color = "black", linewidth = REF_AXIS_TICK_LINEWIDTH),
      axis.ticks.length = unit(REF_AXIS_TICK_LENGTH_PT, "pt"),
      panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.16),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.34),
      legend.position = "none",
      plot.margin = margin(3.0, 8, 1.8, 8)
    )
}

transform_value_ref <- function(x, mode) {
  if (mode == "raw") return(x)
  if (mode == "log10") return(log10(x + PSEUDO))
  stop("Unknown transform mode: ", mode, call. = FALSE)
}

axis_power_ref <- function(metric, transform_mode) {
  if (transform_mode != "raw") return(0)
  if (metric %in% c("cds_snp_density", "cds_deleterious_density")) return(-3)
  0
}

metric_y_label_ref <- function(metric, transform_mode = "raw") {
  base <- switch(
    metric,
    cds_conserved_pct = "Constrained CDS sites (%)",
    cds_snp_density = "CDS SNP density",
    cds_deleterious_density = "Deleterious SNP density",
    cds_deleterious_ratio = "Deleterious SNP ratio",
    metric
  )
  if (transform_mode == "log10") return(bquote(log[10](.(base) + 10^-6)))
  power <- axis_power_ref(metric, transform_mode)
  if (power != 0) return(bquote(.(base)~"(×"~10^.(power)~")"))
  base
}

trim_zeros_ref <- function(x) {
  x <- sub("\\.0+$", "", x)
  sub("(\\.[0-9]*?)0+$", "\\1", x)
}

format_axis_ref <- function(x, metric, transform_mode) {
  if (transform_mode == "log10") {
    return(trim_zeros_ref(format(round(x, 2), trim = TRUE, scientific = FALSE)))
  }
  power <- axis_power_ref(metric, transform_mode)
  z <- as.numeric(x) / (10 ^ power)
  vapply(z, function(v) {
    if (!is.finite(v)) return("")
    if (abs(v) < 1e-12) v <- 0
    av <- abs(v)
    if (abs(v - round(v)) < 1e-8) return(format(round(v), trim = TRUE, scientific = FALSE))
    if (av >= 100) return(trim_zeros_ref(format(round(v, 0), trim = TRUE, scientific = FALSE)))
    if (av >= 10) return(trim_zeros_ref(format(round(v, 1), trim = TRUE, scientific = FALSE, nsmall = 1)))
    if (av >= 1) return(trim_zeros_ref(format(round(v, 2), trim = TRUE, scientific = FALSE, nsmall = 2)))
    if (av >= 0.1) return(trim_zeros_ref(format(round(v, 2), trim = TRUE, scientific = FALSE, nsmall = 2)))
    if (av >= 0.01) return(trim_zeros_ref(format(round(v, 3), trim = TRUE, scientific = FALSE, nsmall = 3)))
    trim_zeros_ref(format(signif(v, 2), trim = TRUE, scientific = FALSE))
  }, character(1))
}

nice_axis_breaks_ref <- function(lims, metric, transform_mode, n = REF_N_Y_BREAKS) {
  vals <- as.numeric(lims)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(seq(0, 1, length.out = n))
  lo <- min(vals)
  hi <- max(vals)
  if (hi <= lo) {
    pad <- max(abs(hi) * 0.05, ifelse(transform_mode == "raw", 1e-8, 0.1))
    lo <- lo - pad
    hi <- hi + pad
  }

  power <- axis_power_ref(metric, transform_mode)
  sf <- 10 ^ power
  lo_s <- lo / sf
  hi_s <- hi / sf

  br_s <- pretty(c(lo_s, hi_s), n = n)
  br_s <- br_s[br_s >= lo_s - 1e-9 & br_s <= hi_s + 1e-9]
  if (length(br_s) != n) br_s <- seq(lo_s, hi_s, length.out = n)
  unique(br_s * sf)
}

get_display_probs_ref <- function(metric, transform_mode) {
  qq <- if (transform_mode == "raw") REF_DISPLAY_QUANTILES_RAW[[metric]] else REF_DISPLAY_QUANTILES_LOG10[[metric]]
  if (is.null(qq) || length(qq) != 2) qq <- c(0, 0.99)
  pmax(0, pmin(1, as.numeric(qq)))
}

get_core_geom_config <- function(geom_mode) {
  if (geom_mode == "box") {
    geom_width <- REF_BOX_WIDTH
    inner_width <- REF_BOX_WIDTH
  } else if (geom_mode == "violin_box") {
    geom_width <- REF_VIOLIN_WIDTH
    inner_width <- REF_VIOLIN_BOX_WIDTH
  } else {
    stop("Unknown geom_mode: ", geom_mode, call. = FALSE)
  }
  list(
    geom_width = geom_width,
    inner_width = inner_width,
    center_step = 2 * REF_BOX_WIDTH
  )
}

get_species_centers_ref <- function(species_in_panel, geom_mode) {
  species_in_panel <- as.character(species_in_panel)
  cfg <- get_core_geom_config(geom_mode)
  centers <- 1 + (seq_along(species_in_panel) - 1) * cfg$center_step
  setNames(centers, species_in_panel)
}

box_whisker_stats_ref <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(lower = NA_real_, q1 = NA_real_, med = NA_real_, q3 = NA_real_, upper = NA_real_))
  }
  qs <- as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE))
  bp <- boxplot.stats(x, coef = 1.5, do.conf = FALSE, do.out = FALSE)$stats
  c(lower = bp[1], q1 = qs[1], med = qs[2], q3 = qs[3], upper = bp[5])
}

compute_species_whisker_df_ref <- function(long_df, metric, transform_mode) {
  long_df %>%
    filter(is.finite(.data$value)) %>%
    mutate(plot_value = transform_value_ref(.data$value, transform_mode)) %>%
    group_by(.data$species) %>%
    summarise(
      lower = box_whisker_stats_ref(.data$plot_value)["lower"],
      q1 = box_whisker_stats_ref(.data$plot_value)["q1"],
      median = box_whisker_stats_ref(.data$plot_value)["med"],
      q3 = box_whisker_stats_ref(.data$plot_value)["q3"],
      upper = box_whisker_stats_ref(.data$plot_value)["upper"],
      .groups = "drop"
    )
}

compute_violin_clip_df_ref <- function(long_df, metric, transform_mode) {
  probs <- get_display_probs_ref(metric, transform_mode)
  long_df %>%
    filter(is.finite(.data$value)) %>%
    mutate(plot_value = transform_value_ref(.data$value, transform_mode)) %>%
    group_by(.data$species) %>%
    summarise(
      violin_lower = safe_quantile_ref(.data$plot_value, probs[1]),
      violin_upper = safe_quantile_ref(.data$plot_value, probs[2]),
      .groups = "drop"
    )
}

build_pairwise_significance_ref <- function(stat_tbl, species_in_panel,
                                            geom_mode, y_base, y_step) {
  if (!SHOW_SIGNIFICANCE || !("group1" %in% names(stat_tbl)) || !("group2" %in% names(stat_tbl))) return(tibble())
  centers <- get_species_centers_ref(species_in_panel, geom_mode)

  out <- stat_tbl %>%
    filter(
      !is.na(.data$group1), !is.na(.data$group2),
      .data$group1 %in% species_in_panel,
      .data$group2 %in% species_in_panel
    ) %>%
    left_join(pairwise_display_order_ref, by = c("group1", "group2"))

  if ("P.adj" %in% names(out)) {
    out <- out %>% mutate(p_use = .data$P.adj)
  } else if ("P.raw" %in% names(out)) {
    out <- out %>% mutate(p_use = .data$P.raw)
  } else {
    out <- out %>% mutate(p_use = NA_real_)
  }

  if ("label" %in% names(out)) {
    out <- out %>% mutate(label = as.character(.data$label))
  } else {
    out <- out %>% mutate(label = vapply(.data$p_use, p_to_star, character(1)))
  }

  out <- out %>%
    mutate(
      xmin = pmin(unname(centers[.data$group1]), unname(centers[.data$group2])),
      xmax = pmax(unname(centers[.data$group1]), unname(centers[.data$group2]))
    ) %>%
    filter(is.finite(.data$p_use), !is.na(.data$ord)) %>%
    distinct(.data$group1, .data$group2, .keep_all = TRUE)

  if (!SHOW_NS) out <- out %>% filter(.data$label != "ns")
  if (nrow(out) == 0) return(tibble())

  out %>%
    arrange(.data$ord) %>%
    mutate(
      y = y_base + (row_number() - 1) * y_step,
      tip = 0.22 * y_step
    ) %>%
    select(.data$xmin, .data$xmax, .data$y, .data$tip, .data$label)
}

add_significance_layers_ref <- function(p, sig) {
  if (nrow(sig) == 0) return(p)
  p +
    geom_segment(
      data = sig,
      aes(x = .data$xmin, xend = .data$xmax, y = .data$y, yend = .data$y),
      inherit.aes = FALSE, linewidth = 0.22, color = "black"
    ) +
    geom_segment(
      data = sig,
      aes(x = .data$xmin, xend = .data$xmin, y = .data$y, yend = .data$y - .data$tip),
      inherit.aes = FALSE, linewidth = 0.22, color = "black"
    ) +
    geom_segment(
      data = sig,
      aes(x = .data$xmax, xend = .data$xmax, y = .data$y, yend = .data$y - .data$tip),
      inherit.aes = FALSE, linewidth = 0.22, color = "black"
    ) +
    geom_text(
      data = sig,
      aes(x = (.data$xmin + .data$xmax) / 2, y = .data$y + 0.18 * .data$tip, label = .data$label),
      inherit.aes = FALSE, size = 2.05, color = "black", vjust = 0
    )
}

compute_shared_axis_ref <- function(long_df, metric, transform_mode,
                                    stat_df = NULL) {
  whisk <- compute_species_whisker_df_ref(long_df, metric, transform_mode)
  if (nrow(whisk) == 0) stop("No finite values available for axis calculation.", call. = FALSE)

  data_lower <- min(whisk$lower, na.rm = TRUE)
  data_upper <- max(whisk$upper, na.rm = TRUE)
  if (!is.finite(data_lower) || !is.finite(data_upper)) stop("Axis whiskers are not finite.", call. = FALSE)
  if (data_upper <= data_lower) {
    pad0 <- max(abs(data_upper) * 0.08, ifelse(transform_mode == "raw", 0.001, 0.1))
    data_lower <- data_lower - pad0
    data_upper <- data_upper + pad0
  }
  data_range <- data_upper - data_lower

  n_brackets <- 3
  if (!is.null(stat_df) && nrow(stat_df) > 0 && all(c("group1", "group2") %in% names(stat_df))) {
    tmp <- stat_df
    if ("label" %in% names(tmp)) {
      tmp <- tmp %>% mutate(lbl = as.character(.data$label))
    } else if ("P.adj" %in% names(tmp)) {
      tmp <- tmp %>% mutate(lbl = vapply(.data$P.adj, p_to_star, character(1)))
    } else if ("P.raw" %in% names(tmp)) {
      tmp <- tmp %>% mutate(lbl = vapply(.data$P.raw, p_to_star, character(1)))
    } else {
      tmp <- tmp %>% mutate(lbl = NA_character_)
    }
    tmp <- tmp %>% left_join(pairwise_display_order_ref, by = c("group1", "group2")) %>% filter(!is.na(.data$ord))
    if (!SHOW_NS) tmp <- tmp %>% filter(.data$lbl != "ns")
    n_brackets <- max(nrow(tmp), 0)
  }

  n_pad <- 0.10 * data_range
  sig_base <- data_upper + 0.05 * data_range
  sig_step <- 0.055 * data_range
  y_lower <- data_lower - n_pad
  y_upper <- if (n_brackets > 0) {
    sig_base + (n_brackets - 1) * sig_step + 1.4 * 0.22 * sig_step + 0.18 * sig_step
  } else {
    data_upper + 0.08 * data_range
  }

  breaks <- nice_axis_breaks_ref(c(data_lower, data_upper), metric, transform_mode, REF_N_Y_BREAKS)
  breaks <- breaks[breaks >= data_lower - 1e-12 & breaks <= data_upper + 1e-12]
  if (length(breaks) != REF_N_Y_BREAKS) breaks <- seq(data_lower, data_upper, length.out = REF_N_Y_BREAKS)

  list(
    data_lower = data_lower,
    data_upper = data_upper,
    data_range = data_range,
    y_limits = c(y_lower, y_upper),
    y_breaks = breaks,
    y_labels = format_axis_ref(breaks, metric, transform_mode),
    sig_base = sig_base,
    sig_step = sig_step,
    n_y = data_lower - 0.065 * data_range
  )
}

plot_core_panel_ref <- function(long_df, stat_tbl, title, metric,
                                geom_mode = "box", transform_mode = "raw",
                                species_in_panel = species_order,
                                scale_mode = "global",
                                show_significance = TRUE,
                                show_n = SHOW_N,
                                y_title = NULL,
                                base_size = 7.4,
                                axis_cfg = NULL) {
  species_in_panel <- as.character(species_in_panel)
  cfg <- get_core_geom_config(geom_mode)
  center_map <- get_species_centers_ref(species_in_panel, geom_mode)

  box_df <- long_df %>%
    filter(as.character(.data$species) %in% species_in_panel, is.finite(.data$value)) %>%
    mutate(
      species_chr = as.character(.data$species),
      xpos = unname(center_map[.data$species_chr]),
      plot_value = transform_value_ref(.data$value, transform_mode)
    ) %>%
    filter(is.finite(.data$plot_value))

  if (nrow(box_df) == 0) stop("No finite plotting values for metric=", metric, ", transform=", transform_mode, call. = FALSE)
  if (anyDuplicated(box_df[, c("SG", "species_chr")])) stop("Plot data contain duplicated SG-species rows for ", metric, ".", call. = FALSE)

  if (is.null(axis_cfg)) {
    axis_cfg <- compute_shared_axis_ref(long_df, metric, transform_mode, stat_tbl)
  }

  violin_df <- box_df
  if (geom_mode == "violin_box") {
    clip_df <- compute_violin_clip_df_ref(long_df, metric, transform_mode) %>% rename(species_chr = .data$species)
    violin_df <- violin_df %>%
      left_join(clip_df, by = "species_chr") %>%
      mutate(
        violin_lower = ifelse(is.finite(.data$violin_lower), .data$violin_lower, axis_cfg$data_lower),
        violin_upper = ifelse(is.finite(.data$violin_upper), .data$violin_upper, axis_cfg$data_upper),
        plot_value_violin = pmin(pmax(.data$plot_value, .data$violin_lower), .data$violin_upper)
      )
  }

  sig <- if (show_significance) {
    build_pairwise_significance_ref(stat_tbl, species_in_panel, geom_mode, axis_cfg$sig_base, axis_cfg$sig_step)
  } else tibble()

  n_df <- box_df %>%
    count(.data$species_chr, .data$xpos, name = "n") %>%
    mutate(y = axis_cfg$n_y, label = paste0("n = ", .data$n, " SGs"))

  x_breaks <- unname(center_map)
  x_labels <- parse(text = unname(species_axis_display[species_in_panel]))
  x_pad <- cfg$geom_width * 1.6
  x_min <- min(x_breaks) - x_pad
  x_max <- max(x_breaks) + x_pad

  p <- ggplot()
  if (geom_mode == "violin_box") {
    p <- p +
      geom_violin(
        data = violin_df,
        aes(x = .data$xpos, y = .data$plot_value_violin, fill = .data$species, color = .data$species, group = .data$species),
        width = cfg$geom_width, scale = "width", trim = TRUE,
        linewidth = 0.22, alpha = 0.90
      ) +
      geom_boxplot(
        data = box_df,
        aes(x = .data$xpos, y = .data$plot_value, fill = .data$species, group = .data$species),
        width = cfg$inner_width, coef = 1.5,
        fill = "white", color = "black", outlier.shape = NA,
        linewidth = 0.30
      )
  } else {
    p <- p +
      geom_boxplot(
        data = box_df,
        aes(x = .data$xpos, y = .data$plot_value, fill = .data$species, color = .data$species, group = .data$species),
        width = cfg$geom_width, coef = 1.5,
        outlier.shape = NA, linewidth = 0.28, alpha = 0.88
      )
  }

  p <- p +
    stat_summary(
      data = box_df,
      aes(x = .data$xpos, y = .data$plot_value, group = .data$species),
      fun = median, geom = "point", shape = 21,
      size = 1.10, stroke = 0.24,
      fill = "black", color = "black"
    )

  p <- add_significance_layers_ref(p, sig)

  if (show_n) {
    p <- p +
      geom_text(
        data = n_df,
        aes(x = .data$xpos, y = .data$y, label = .data$label),
        inherit.aes = FALSE, size = 1.70, hjust = 0.5, vjust = 0, color = "#333333"
      )
  }

  if (is.null(y_title)) y_title <- metric_y_label_ref(metric, transform_mode)

  p +
    scale_fill_manual(values = species_fill, limits = species_order, guide = "none") +
    scale_color_manual(values = species_line, limits = species_order, guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = axis_cfg$y_breaks, labels = axis_cfg$y_labels, expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = axis_cfg$y_limits, clip = "off") +
    labs(x = NULL, y = y_title, title = title) +
    theme_core_reference(base_size = base_size)
}

combine_row_shared_y_ref <- function(plot_list, y_title, tag_levels = NULL) {
  if (length(plot_list) == 0) stop("No plots supplied.", call. = FALSE)
  plot_list <- lapply(seq_along(plot_list), function(i) {
    if (i == 1) plot_list[[i]] + labs(y = y_title)
    else plot_list[[i]] + labs(y = NULL) + theme(axis.title.y = element_blank())
  })
  out <- wrap_plots(plot_list, nrow = 1)
  if (!is.null(tag_levels)) out <- out + plot_annotation(tag_levels = tag_levels)
  out
}

plot_threshold_row_ref <- function(long_all, stats_all, metric, data_type,
                                   geom_mode, transform_mode,
                                   scale_mode = "global") {
  row_long <- long_all %>%
    filter(.data$metric == .env$metric, .data$data_type == .env$data_type)
  row_stats <- stats_all %>%
    filter(.data$metric == .env$metric, .data$data_type == .env$data_type)
  axis_cfg <- compute_shared_axis_ref(row_long, metric, transform_mode, row_stats)

  plots <- lapply(THRESHOLDS, function(thr) {
    long_sub <- row_long %>% filter(as.character(.data$threshold) == .env$thr)
    stat_sub <- row_stats %>% filter(as.character(.data$threshold) == .env$thr)
    plot_core_panel_ref(
      long_sub, stat_sub,
      title = paste0("GERP > ", sub("gt", "", thr)),
      metric = metric,
      geom_mode = geom_mode,
      transform_mode = transform_mode,
      species_in_panel = species_order,
      scale_mode = scale_mode,
      show_significance = TRUE,
      show_n = SHOW_N,
      y_title = NULL,
      axis_cfg = axis_cfg
    )
  })
  combine_row_shared_y_ref(plots, metric_y_label_ref(metric, transform_mode))
}

plot_snp_density_combined_ref <- function(snp_long, snp_stats, data_type,
                                          geom_mode, transform_mode) {
  long_sub <- snp_long %>%
    filter(
      as.character(.data$threshold) == SNP_DENSITY_THRESHOLD,
      .data$metric == "cds_snp_density",
      .data$data_type == .env$data_type
    )
  stat_sub <- snp_stats %>%
    filter(
      as.character(.data$threshold) == SNP_DENSITY_THRESHOLD,
      .data$metric == "cds_snp_density",
      .data$data_type == .env$data_type
    )
  axis_cfg <- compute_shared_axis_ref(long_sub, "cds_snp_density", transform_mode, stat_sub)
  plot_core_panel_ref(
    long_sub, stat_sub,
    title = "CDS SNP density",
    metric = "cds_snp_density",
    geom_mode = geom_mode,
    transform_mode = transform_mode,
    species_in_panel = species_order,
    scale_mode = "global",
    show_significance = TRUE,
    show_n = SHOW_N,
    axis_cfg = axis_cfg
  )
}

plot_snp_density_split_ref <- function(snp_long, data_type,
                                       geom_mode, transform_mode) {
  long_sub <- snp_long %>%
    filter(
      as.character(.data$threshold) == SNP_DENSITY_THRESHOLD,
      .data$metric == "cds_snp_density",
      .data$data_type == .env$data_type
    )
  axis_cfg <- compute_shared_axis_ref(long_sub, "cds_snp_density", transform_mode, NULL)
  plots <- lapply(species_order, function(sp) {
    plot_core_panel_ref(
      long_sub %>% filter(as.character(.data$species) == .env$sp),
      tibble(),
      title = parse(text = species_axis_display[[sp]]),
      metric = "cds_snp_density",
      geom_mode = geom_mode,
      transform_mode = transform_mode,
      species_in_panel = sp,
      scale_mode = "global",
      show_significance = FALSE,
      show_n = SHOW_N,
      y_title = NULL,
      axis_cfg = axis_cfg
    )
  })
  combine_row_shared_y_ref(plots, metric_y_label_ref("cds_snp_density", transform_mode))
}

plot_snp0_panel_ref <- function(summary_df, stat_tbl) {
  cfg <- get_core_geom_config("box")
  center_map <- get_species_centers_ref(species_order, "box")
  df <- summary_df %>% mutate(species_chr = as.character(.data$species), xpos = unname(center_map[.data$species_chr]))

  sig_base <- 1.03
  sig_step <- 0.05
  sig <- build_pairwise_significance_ref(stat_tbl, species_order, "box", sig_base, sig_step)
  y_upper <- if (nrow(sig) > 0) max(sig$y + 1.35 * sig$tip + 0.15 * sig_step) else 1.06

  p <- ggplot(df, aes(x = .data$xpos, y = .data$prop_SNP0, fill = .data$species, color = .data$species)) +
    geom_col(width = REF_SNP0_BAR_WIDTH, linewidth = 0.28, alpha = 0.88) +
    geom_text(aes(label = paste0(.data$n_SNP0_SG, "/", .data$n_complete_SG)), vjust = -0.35, size = 1.72, color = "black")
  p <- add_significance_layers_ref(p, sig)

  x_breaks <- unname(center_map)
  x_pad <- cfg$geom_width * 1.6
  x_min <- min(x_breaks) - x_pad
  x_max <- max(x_breaks) + x_pad

  p +
    scale_fill_manual(values = species_fill, limits = species_order, guide = "none") +
    scale_color_manual(values = species_line, limits = species_order, guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = parse(text = unname(species_axis_display[species_order])), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = seq(0, 1, length.out = REF_N_Y_BREAKS), labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, y_upper), clip = "off") +
    labs(x = NULL, y = "SNP-free strict-core SGs (%)", title = "Total CDS SNP count = 0") +
    theme_core_reference()
}

center_single_panel_row_ref <- function(...) {
  panels <- list(...)
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) == 1) {
    return(wrap_plots(plot_spacer(), panels[[1]], plot_spacer(), nrow = 1, widths = c(1, 1, 1)))
  }
  wrap_plots(panels, nrow = 1)
}

make_a4_core_figure_ref <- function(
    data_type, transform_mode, geom_mode,
    gerp_plot_long, gerp_plot_stats,
    snp_long, snp_stats,
    snp0_panel,
    extra_bottom_panel = NULL) {

  p_gerp <- plot_threshold_row_ref(
    gerp_plot_long, gerp_plot_stats,
    metric = "cds_conserved_pct", data_type = data_type,
    geom_mode = geom_mode, transform_mode = transform_mode
  )

  p_snp <- plot_snp_density_combined_ref(snp_long, snp_stats, data_type, geom_mode, transform_mode)
  p_snp_row <- center_single_panel_row_ref(p_snp)

  p_del <- plot_threshold_row_ref(
    snp_long, snp_stats,
    metric = "cds_deleterious_density", data_type = data_type,
    geom_mode = geom_mode, transform_mode = transform_mode
  )

  p_ratio <- plot_threshold_row_ref(
    snp_long, snp_stats,
    metric = "cds_deleterious_ratio", data_type = "SNPpos",
    geom_mode = geom_mode, transform_mode = transform_mode
  )

  bottom_row <- center_single_panel_row_ref(snp0_panel, extra_bottom_panel)

  if (data_type == "ALL") {
    out <- p_gerp / p_snp_row / p_del / p_ratio / bottom_row +
      plot_layout(heights = c(1, 1, 1, 1, 1))
  } else {
    out <- p_gerp / p_snp_row / p_del / p_ratio +
      plot_layout(heights = c(1, 1, 1, 1))
  }

  out +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 10, face = "bold"),
        plot.margin = margin(3.5, 8, 3.5, 8)
      )
    )
}


## ============================================================================
## 9. Main workflow

## ============================================================================

## ============================================================================

check_all_inputs()
pan_map <- read_sg_pan(PAN_FILE)

## SNP and deleterious-SNP tables ---------------------------------------------
snp_wide <- setNames(lapply(THRESHOLDS, function(thr) {
  message("[INFO] Reading SNP/deleterious inputs: ", thr)
  build_snp_wide_one_threshold(thr, pan_map)
}), THRESHOLDS)

snp_input_audit <- audit_snp_inputs_across_thresholds(snp_wide)
write_tsv_base(
  snp_input_audit,
  file.path(
    TABLE_DIR,
    "SNP_total_and_CDS_length_threshold_consistency_audit.tsv"
  )
)
if (any(
  snp_input_audit$n_missing_pattern_different > 0 |
    snp_input_audit$n_value_different > 0,
  na.rm = TRUE
)) {
  warning(
    "Total SNP counts, CDS lengths, or missingness differ among ",
    "threshold-specific inputs. Threshold-independent SNP density and SNP0 ",
    "are calculated only from ", SNP_DENSITY_THRESHOLD,
    ". Inspect the consistency audit.",
    call. = FALSE
  )
}

write_tsv_base(
  bind_rows(snp_wide),
  file.path(
    TABLE_DIR, "Strict_core_SG_SNP_deleterious_wide_all_thresholds.tsv"
  )
)
write_tsv_base(
  bind_rows(lapply(snp_wide, ratio_audit)),
  file.path(TABLE_DIR, "Deleterious_ratio_input_audit.tsv")
)

## SNP density is threshold-independent and is therefore calculated once.
snp_density_long <- bind_rows(
  metric_long_from_wide(
    snp_wide[[SNP_DENSITY_THRESHOLD]], "cds_snp_density", "ALL"
  ),
  metric_long_from_wide(
    snp_wide[[SNP_DENSITY_THRESHOLD]], "cds_snp_density", "SNPpos"
  )
)

## Deleterious metrics depend on the GERP cutoff.
deleterious_long <- bind_rows(lapply(THRESHOLDS, function(thr) {
  bind_rows(
    metric_long_from_wide(
      snp_wide[[thr]], "cds_deleterious_density", "ALL"
    ),
    metric_long_from_wide(
      snp_wide[[thr]], "cds_deleterious_density", "SNPpos"
    ),
    metric_long_from_wide(
      snp_wide[[thr]], "cds_deleterious_ratio", "SNPpos"
    )
  )
}))

snp_long_all <- bind_rows(snp_density_long, deleterious_long)
write_tsv_base(
  snp_long_all,
  file.path(TABLE_DIR, "Strict_core_SNP_deleterious_metrics_long.tsv")
)
write_tsv_base(
  summarise_long_metric(snp_long_all),
  file.path(TABLE_DIR, "Strict_core_SNP_deleterious_metric_summary.tsv")
)

snp_density_stats <- bind_rows(
  paired_continuous_stats(
    snp_wide[[SNP_DENSITY_THRESHOLD]], "cds_snp_density",
    SNP_DENSITY_THRESHOLD, "ALL"
  ),
  paired_continuous_stats(
    snp_wide[[SNP_DENSITY_THRESHOLD]], "cds_snp_density",
    SNP_DENSITY_THRESHOLD, "SNPpos"
  )
)

deleterious_stats <- bind_rows(lapply(THRESHOLDS, function(thr) {
  bind_rows(
    paired_continuous_stats(
      snp_wide[[thr]], "cds_deleterious_density", thr, "ALL"
    ),
    paired_continuous_stats(
      snp_wide[[thr]], "cds_deleterious_density", thr, "SNPpos"
    ),
    paired_continuous_stats(
      snp_wide[[thr]], "cds_deleterious_ratio", thr, "SNPpos"
    )
  )
}))

snp_stats_all <- bind_rows(snp_density_stats, deleterious_stats)
write_tsv_base(
  snp_stats_all,
  file.path(STAT_DIR, "Strict_core_SNP_deleterious_paired_stats.tsv")
)

## SNP0 is independent of the GERP cutoff and is calculated once from gt4.
snp0_obj <- snp0_paired_stats(
  snp_wide[[SNP_DENSITY_THRESHOLD]], SNP_DENSITY_THRESHOLD
)
snp0_summary_all <- snp0_obj$summary
snp0_stats_all <- snp0_obj$stats
write_tsv_base(
  snp0_summary_all,
  file.path(TABLE_DIR, "Strict_core_SNP0_SG_proportion.tsv")
)
write_tsv_base(
  snp0_stats_all,
  file.path(STAT_DIR, "Strict_core_SNP0_CochranQ_exact_McNemar.tsv")
)

## GERP constrained-CDS tables ------------------------------------------------
gerp_wide <- setNames(lapply(THRESHOLDS, function(thr) {
  message("[INFO] Reading GERP CDS inputs: ", thr)
  build_gerp_wide_one_threshold(thr, pan_map)
}), THRESHOLDS)

write_tsv_base(
  bind_rows(gerp_wide),
  file.path(
    TABLE_DIR, "Strict_core_SG_GERP_CDS_wide_all_thresholds.tsv"
  )
)

gerp_long_all <- bind_rows(lapply(gerp_wide, gerp_long_from_wide))
write_tsv_base(
  gerp_long_all,
  file.path(TABLE_DIR, "Strict_core_GERP_CDS_constrained_pct_long.tsv")
)
write_tsv_base(
  summarise_long_metric(gerp_long_all),
  file.path(TABLE_DIR, "Strict_core_GERP_CDS_constrained_pct_summary.tsv")
)

gerp_stats_all <- bind_rows(lapply(THRESHOLDS, function(thr) {
  paired_continuous_stats(
    gerp_wide[[thr]], "cds_conserved_pct", thr, "GERP_CDS"
  )
}))
write_tsv_base(
  gerp_stats_all,
  file.path(STAT_DIR, "Strict_core_GERP_CDS_paired_stats.tsv")
)

flush_qc_summaries()

## Exact number of matched SGs used by both plots and tests.
paired_n_summary <- bind_rows(
  snp_long_all %>%
    distinct(.data$threshold, .data$data_type, .data$metric, .data$SG) %>%
    count(
      .data$threshold, .data$data_type, .data$metric,
      name = "n_paired_SG"
    ),
  gerp_long_all %>%
    distinct(.data$threshold, .data$data_type, .data$metric, .data$SG) %>%
    count(
      .data$threshold, .data$data_type, .data$metric,
      name = "n_paired_SG"
    )
)
write_tsv_base(
  paired_n_summary,
  file.path(TABLE_DIR, "Strict_core_metric_paired_sample_sizes.tsv")
)

## Figures --------------------------------------------------------------------

## Canonical SNP-positive strict-SG set. Total CDS SNP counts should be
## threshold-independent; gt4 is retained as the audited canonical source.
canonical_snp_positive_sg <- prepare_paired_metric_wide(
  snp_wide[[SNP_DENSITY_THRESHOLD]],
  metric = "cds_snp_density",
  data_type = "SNPpos"
)$SG

## GERP plotting datasets are emitted for both ALL and SNPpos A4 figures.
## The SNPpos GERP row is restricted to the same strict SGs that are SNP-positive
## in all three species, then complete cases for the GERP metric are retained.
gerp_plot_objects <- lapply(CONTINUOUS_DATA_TYPES, function(dt) {
  long_obj <- bind_rows(lapply(THRESHOLDS, function(thr) {
    wide_use <- gerp_wide[[thr]]
    if (dt == "SNPpos") {
      wide_use <- wide_use %>%
        filter(.data$SG %in% canonical_snp_positive_sg)
    }
    gerp_long_from_wide(wide_use) %>%
      mutate(data_type = dt)
  }))

  stat_obj <- bind_rows(lapply(THRESHOLDS, function(thr) {
    wide_use <- gerp_wide[[thr]]
    if (dt == "SNPpos") {
      wide_use <- wide_use %>%
        filter(.data$SG %in% canonical_snp_positive_sg)
    }
    paired_continuous_stats(
      wide_use, "cds_conserved_pct", thr, "GERP_CDS"
    ) %>%
      mutate(data_type = dt)
  }))

  list(long = long_obj, stats = stat_obj)
})
names(gerp_plot_objects) <- CONTINUOUS_DATA_TYPES

gerp_plot_long <- bind_rows(lapply(gerp_plot_objects, function(x) x$long))
gerp_plot_stats <- bind_rows(lapply(gerp_plot_objects, function(x) x$stats))

write_tsv_base(
  gerp_plot_long,
  file.path(TABLE_DIR, "Strict_core_GERP_CDS_plot_data_ALL_SNPpos.tsv")
)
write_tsv_base(
  gerp_plot_stats,
  file.path(STAT_DIR, "Strict_core_GERP_CDS_plot_stats_ALL_SNPpos.tsv")
)
write_tsv_base(
  tibble(
    SNP_status_definition = "SNP-positive in all three species using canonical gt4 total-CDS-SNP input",
    n_canonical_SNP_positive_strict_SG = length(canonical_snp_positive_sg)
  ),
  file.path(TABLE_DIR, "Strict_core_SNPpos_subset_definition.tsv")
)

p_snp0_base <- plot_snp0_panel_ref(snp0_summary_all, snp0_stats_all)
safe_ggsave(
  file.path(FIG_DIR, "Fig_strict_core_SNP0_proportion_reference_style.pdf"),
  p_snp0_base, ROW_WIDTH / 3, ROW_HEIGHT, FIG_DPI
)
safe_ggsave(
  file.path(FIG_DIR, "Fig_strict_core_SNP0_proportion_reference_style.png"),
  p_snp0_base, ROW_WIDTH / 3, ROW_HEIGHT, FIG_DPI
)

for (geom_mode in GEOM_MODES_TO_SAVE) {
  for (transform_mode in DENSITY_TRANSFORMS) {
    mode_tag <- ifelse(transform_mode == "raw", "raw_whiskerY", "log10_whiskerY")

    ## GERP, SNP density, deleterious density, and full A4 figures are written
    ## separately for ALL and for the all-three-species SNP-positive subset.
    for (dt in CONTINUOUS_DATA_TYPES) {
      p_gerp_row <- plot_threshold_row_ref(
        gerp_plot_long, gerp_plot_stats,
        metric = "cds_conserved_pct", data_type = dt,
        geom_mode = geom_mode, transform_mode = transform_mode
      )
      prefix <- paste0(
        "Fig_strict_core_GERP_CDS_", dt, "_gt2_gt4_gt6_",
        mode_tag, "_", geom_mode
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".pdf")),
        p_gerp_row + plot_annotation(tag_levels = "A"),
        ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".png")),
        p_gerp_row + plot_annotation(tag_levels = "A"),
        ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )

      ## Three species in one combined panel with species-local Y axes.
      p_snp_combined <- plot_snp_density_combined_ref(
        snp_long_all, snp_stats_all, dt,
        geom_mode, transform_mode
      )
      prefix <- paste0(
        "Fig_strict_core_SNP_density_combined_species_", dt, "_",
        mode_tag, "_", geom_mode
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".pdf")),
        p_snp_combined, ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".png")),
        p_snp_combined, ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )

      ## Three species shown as three separate panels. This is useful for
      ## inspecting each distribution on a conventional single Y axis.
      p_snp_split <- plot_snp_density_split_ref(
        snp_long_all, dt, geom_mode, transform_mode
      )
      prefix <- paste0(
        "Fig_strict_core_SNP_density_split_three_species_", dt, "_",
        mode_tag, "_", geom_mode
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".pdf")),
        p_snp_split, ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".png")),
        p_snp_split, ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )

      p_del_row <- plot_threshold_row_ref(
        snp_long_all, snp_stats_all,
        metric = "cds_deleterious_density", data_type = dt,
        geom_mode = geom_mode, transform_mode = transform_mode
      )
      prefix <- paste0(
        "Fig_strict_core_deleterious_density_", dt,
        "_gt2_gt4_gt6_", mode_tag, "_", geom_mode
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".pdf")),
        p_del_row + plot_annotation(tag_levels = "A"),
        ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".png")),
        p_del_row + plot_annotation(tag_levels = "A"),
        ROW_WIDTH, ROW_HEIGHT, FIG_DPI
      )

      p_a4 <- make_a4_core_figure_ref(
        data_type = dt,
        transform_mode = transform_mode,
        geom_mode = geom_mode,
        gerp_plot_long = gerp_plot_long,
        gerp_plot_stats = gerp_plot_stats,
        snp_long = snp_long_all,
        snp_stats = snp_stats_all,
        snp0_panel = p_snp0_base
      )
      prefix <- paste0(
        "Fig_strict_core_3species_A4_", dt, "_",
        mode_tag, "_", geom_mode
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".pdf")),
        p_a4, A4_WIDTH, A4_HEIGHT, FIG_DPI
      )
      safe_ggsave(
        file.path(FIG_DIR, paste0(prefix, ".png")),
        p_a4, A4_WIDTH, A4_HEIGHT, FIG_DPI
      )
    }

    ## Ratio is defined only for SNP-positive matched SGs; it is output once
    ## per transform/geometry rather than duplicated under an ALL label.
    p_ratio_row <- plot_threshold_row_ref(
      snp_long_all, snp_stats_all,
      metric = "cds_deleterious_ratio", data_type = "SNPpos",
      geom_mode = geom_mode, transform_mode = transform_mode
    )
    prefix <- paste0(
      "Fig_strict_core_deleterious_ratio_SNPpos_gt2_gt4_gt6_",
      mode_tag, "_", geom_mode
    )
    safe_ggsave(
      file.path(FIG_DIR, paste0(prefix, ".pdf")),
      p_ratio_row + plot_annotation(tag_levels = "A"),
      ROW_WIDTH, ROW_HEIGHT, FIG_DPI
    )
    safe_ggsave(
      file.path(FIG_DIR, paste0(prefix, ".png")),
      p_ratio_row + plot_annotation(tag_levels = "A"),
      ROW_WIDTH, ROW_HEIGHT, FIG_DPI
    )
  }
}

message("[DONE] Shared-Y whisker-based strict-core figures written to: ", FIG_DIR)
message("[DONE] A4 combinations: 2 data types × 2 transforms × 2 geometries.")
message("[DONE] Results written to: ", OUTDIR)
message("[DONE] Input QC tables written to: ", QC_DIR)
