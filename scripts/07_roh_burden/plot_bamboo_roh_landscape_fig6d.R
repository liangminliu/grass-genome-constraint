#!/usr/bin/env Rscript
# Fig. 6D: moso bamboo only. Adapted from the executed bamboo ROH landscape.
# Usage: Rscript plot_bamboo_roh_landscape_fig6d.R <S3.xlsx> <Ped_ROH_FROH_table_prefix> <output_prefix>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript plot_bamboo_roh_landscape_fig6d.R <S3.xlsx> <Ped_ROH_FROH_table_prefix> <output_prefix>")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install R package 'readxl'")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install R package 'ggplot2'")

sheet <- as.data.frame(readxl::read_excel(args[[1]], sheet = "S3B_bamboo_Ho_ROH", skip = 3))
sheet <- sheet[!is.na(sheet[[3]]) & trimws(as.character(sheet[[3]])) != "", , drop = FALSE]
ids <- as.character(sheet[[3]])
if (length(ids) != 63L || anyDuplicated(ids)) stop("Expected 63 unique final bamboo samples in S3B")

prefix <- args[[2]]
read_input <- function(suffix) {
  path <- paste0(prefix, suffix)
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
froh <- read_input(".froh.tsv")
chr <- read_input(".chr_info.tsv")
roh <- read_input(".roh_segments.clean.tsv")
stopifnot(all(c("sample", "FROH") %in% names(froh)),
          all(c("chr", "chr_start", "chr_end", "chr_len", "offset", "chr_index") %in% names(chr)),
          all(c("sample", "chr", "start_bp", "end_bp", "length_bp") %in% names(roh)))
froh$sample <- as.character(froh$sample)
roh$sample <- as.character(roh$sample)
chr$chr <- as.character(chr$chr)
roh$chr <- as.character(roh$chr)
froh <- froh[froh$sample %in% ids, , drop = FALSE]
roh <- roh[roh$sample %in% ids & roh$length_bp >= 100000, , drop = FALSE]
if (nrow(froh) != 63L || anyDuplicated(froh$sample) || !setequal(froh$sample, ids)) {
  stop("Bamboo FROH table does not match the 63 approved S3B samples")
}
if (!all(is.finite(froh$FROH)) || any(froh$FROH < 0 | froh$FROH > 1)) stop("Invalid FROH values")
if (!all(roh$chr %in% chr$chr)) stop("ROH chromosome absent from chromosome information")

# The executed bamboo plot displayed the five highest and five lowest non-PE
# samples. PE-prefixed samples remain in Fig. 6A-C and in all summary metrics.
eligible <- froh[!grepl("^PE", froh$sample), , drop = FALSE]
eligible <- eligible[order(eligible$FROH, eligible$sample), , drop = FALSE]
if (nrow(eligible) < 10L) stop("Fewer than ten eligible display samples")
selected <- rbind(eligible[seq.int(nrow(eligible), nrow(eligible) - 4L), , drop = FALSE],
                  eligible[seq_len(5L), , drop = FALSE])
selected$track <- rev(seq_len(nrow(selected)))
selected$label <- sprintf("%s (%.4f)", selected$sample, selected$FROH)

roh <- roh[roh$sample %in% selected$sample, , drop = FALSE]
if (!nrow(roh)) stop("No ROH segments for selected samples")
ix <- match(roh$chr, chr$chr)
roh$track <- selected$track[match(roh$sample, selected$sample)]
start <- pmax(roh$start_bp, chr$chr_start[ix])
end <- pmin(roh$end_bp, chr$chr_end[ix])
if (any(end < start)) stop("ROH coordinates fall outside chromosome bounds")
roh$xstart <- chr$offset[ix] + start - chr$chr_start[ix]
roh$xend <- chr$offset[ix] + end - chr$chr_start[ix]
mid <- (roh$xstart + roh$xend) / 2
# Minimum width is a visual aid only; unexpanded coordinates remain in the data.
visual_width <- pmax(roh$xend - roh$xstart, 1500000)
roh$xmin <- pmax(chr$offset[ix], mid - visual_width / 2)
roh$xmax <- pmin(chr$offset[ix] + chr$chr_len[ix], mid + visual_width / 2)
roh$shade <- ifelse(chr$chr_index[ix] %% 2L == 1L, "odd", "even")
chr$xmin <- chr$offset
chr$xmax <- chr$offset + chr$chr_len
chr$center <- (chr$xmin + chr$xmax) / 2
chr$shade <- ifelse(chr$chr_index %% 2L == 1L, "odd", "even")

baseline <- data.frame(track = selected$track,
                       xmin = min(chr$xmin), xmax = max(chr$xmax))
p <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = chr, ggplot2::aes(xmin = xmin, xmax = xmax,
                      ymin = -Inf, ymax = Inf, fill = shade), alpha = 0.07) +
  ggplot2::geom_rect(data = baseline, ggplot2::aes(xmin = xmin, xmax = xmax,
                      ymin = track - 0.045, ymax = track + 0.045), fill = "grey65") +
  ggplot2::geom_rect(data = roh, ggplot2::aes(xmin = xmin, xmax = xmax,
                      ymin = track - 0.24, ymax = track + 0.24, fill = shade), alpha = 0.92) +
  ggplot2::scale_fill_manual(values = c(odd = "#2C5C65", even = "#79A8A7"), guide = "none") +
  ggplot2::scale_x_continuous(breaks = chr$center, labels = chr$chr_label,
                              expand = ggplot2::expansion(mult = 0.005)) +
  ggplot2::scale_y_continuous(breaks = selected$track, labels = selected$label,
                              expand = ggplot2::expansion(mult = 0.04)) +
  ggplot2::labs(x = "Chromosome", y = "Bamboo individual (FROH)") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
                 axis.text.y = ggplot2::element_text(size = 7))

out <- args[[3]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(paste0(out, ".pdf"), width = 13.5, height = 5.2, family = "sans")
print(p)
grDevices::dev.off()
grDevices::png(paste0(out, ".png"), width = 8100, height = 3120, res = 600, type = "cairo")
print(p)
grDevices::dev.off()
write.table(selected[, c("sample", "FROH")], paste0(out, ".displayed_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote bamboo Fig. 6D PDF, PNG and selected-sample TSV")
