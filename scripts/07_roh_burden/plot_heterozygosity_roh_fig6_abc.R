#!/usr/bin/env Rscript
# Standalone recreation of Fig. 6A-C from the final S3B-D sample tables.
# Usage: Rscript plot_heterozygosity_roh_fig6_abc.R <S3.xlsx> <output_prefix>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_heterozygosity_roh_fig6_abc.R <S3.xlsx> <output_prefix>")
}
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the R package 'readxl'")

input <- args[[1]]
prefix <- args[[2]]
if (!file.exists(input)) stop("Input workbook not found: ", input)
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

panels <- list(
  list(sheet = "S3B_bamboo_Ho_ROH", label = "Moso bamboo", n = 63L,
       class_columns = c(15L, 18L, 21L, 24L)),
  list(sheet = "S3C_teosinte_Ho_ROH", label = "Teosinte", n = 90L,
       class_columns = c(16L, 19L, 22L, 25L)),
  list(sheet = "S3D_Aegilops_Ho_ROH", label = "Tausch's goatgrass", n = 140L,
       class_columns = c(15L, 18L, 21L, 24L))
)

read_panel <- function(spec) {
  x <- as.data.frame(readxl::read_excel(input, sheet = spec$sheet, skip = 3))
  x <- x[!is.na(x[[3]]) & trimws(as.character(x[[3]])) != "", , drop = FALSE]
  if (nrow(x) != spec$n || anyDuplicated(x[[3]])) {
    stop("Unexpected sample count or duplicate ID in ", spec$sheet)
  }
  ho <- as.numeric(x[[9]])
  froh <- as.numeric(x[[10]])
  classes <- vapply(spec$class_columns, function(j) as.numeric(x[[j]]), numeric(nrow(x)))
  if (any(!is.finite(c(ho, froh, classes)))) stop("Missing/non-numeric values in ", spec$sheet)
  if (any(ho < 0 | ho > 1) || any(froh < 0 | froh > 1) || any(classes < 0)) {
    stop("Out-of-range values in ", spec$sheet)
  }
  if (any(abs(rowSums(classes) - froh) > 1e-7)) {
    stop("ROH-class fractions do not sum to FROH in ", spec$sheet)
  }
  list(ho = ho, froh = froh, cumulative = t(apply(classes, 1, cumsum)),
       label = paste0(spec$label, " (n = ", spec$n, ")"))
}
data <- lapply(panels, read_panel)

draw <- function() {
  layout(matrix(1:9, nrow = 3, byrow = TRUE), heights = c(0.85, 0.85, 1.25))
  par(family = "sans", mar = c(3.1, 3.5, 2.1, 0.9), mgp = c(2.0, 0.55, 0),
      tcl = -0.25, cex = 0.76, cex.axis = 0.76, cex.lab = 0.78)
  for (i in seq_along(data)) {
    d <- data[[i]]
    hist(d$ho, breaks = pretty(range(d$ho), n = 18), col = "#ACB5D2",
         border = "grey65", main = d$label,
         xlab = expression("Observed heterozygosity (" * italic(H)[o] * ")"),
         ylab = if (i == 1L) "Individuals" else "", cex.main = 0.84)
    if (i == 1L) mtext("A", side = 3, adj = -0.12, line = 0.35, font = 2, cex = 1.2)
  }
  for (i in seq_along(data)) {
    d <- data[[i]]
    hist(d$froh, breaks = pretty(range(d$froh), n = 18), col = "#BDC2B6",
         border = "grey65", main = "", xlab = expression(italic(F)[ROH]),
         ylab = if (i == 1L) "Individuals" else "")
    if (i == 1L) mtext("B", side = 3, adj = -0.12, line = 0.35, font = 2, cex = 1.2)
  }
  par(mar = c(4.2, 3.5, 1.5, 0.9))
  for (i in seq_along(data)) {
    z <- data[[i]]$cumulative
    plot(NA, xlim = c(0.8, 4.2), ylim = c(0, max(z) * 1.08), xaxt = "n",
         xlab = "ROH length class",
         ylab = if (i == 1L) expression("Cumulative " * italic(F)[ROH]) else "",
         bty = "l")
    axis(1, at = 1:4, labels = c("100–250 kb", "250–500 kb", "500 kb–1 Mb", ">1 Mb"),
         cex.axis = 0.68)
    for (j in seq_len(nrow(z))) lines(1:4, z[j, ], col = grDevices::adjustcolor("grey30", 0.23), lwd = 0.65)
    lines(1:4, colMeans(z), col = "black", lwd = 2)
    points(1:4, colMeans(z), col = "black", pch = 16, cex = 0.7)
    if (i == 1L) mtext("C", side = 3, adj = -0.12, line = 0.35, font = 2, cex = 1.2)
  }
}

grDevices::pdf(paste0(prefix, ".pdf"), width = 17.9, height = 11.4, family = "sans")
draw()
grDevices::dev.off()
grDevices::png(paste0(prefix, ".png"), width = 7160, height = 4560, res = 400, type = "cairo")
draw()
grDevices::dev.off()
message("Wrote ", prefix, ".pdf and .png")
