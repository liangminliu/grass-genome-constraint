#!/usr/bin/env Rscript
# Fig. 3A/B and Fig. S1: chromosome, GERP, gene, TE, SNP and pi tracks.
# Usage: Rscript plot_circular_tracks.R <500kb_tracks.tsv> <output.pdf> <species_label>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript plot_circular_tracks.R <500kb_tracks.tsv> <output.pdf> <species_label>")
if (!requireNamespace("circlize", quietly = TRUE)) stop("Install R package 'circlize'")
d <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
need <- c("chrom", "start", "end", "gerp_mean", "gene_fraction", "te_fraction", "snp_per_kb", "pi_mean")
if (!all(need %in% names(d))) stop("Missing required 500-kb track columns")
if (!nrow(d) || any(d$end <= d$start) || anyDuplicated(d[c("chrom", "start")])) stop("Invalid window table")
chrom <- unique(d$chrom)
lengths <- vapply(chrom, function(x) max(d$end[d$chrom == x]), numeric(1))
tracks <- c("gerp_mean", "gene_fraction", "te_fraction", "snp_per_kb", "pi_mean")
colors <- c("#304A69", "#527C8A", "#A68569", "#924E4B", "#59754D")
names(colors) <- tracks

rescale <- function(x, name) {
  if (name %in% c("gene_fraction", "te_fraction")) return(pmin(pmax(x, 0), 1))
  finite <- x[is.finite(x)]
  if (!length(finite)) return(rep(NA_real_, length(x)))
  lo <- min(finite)
  hi <- max(finite)
  if (hi == lo) return(rep(0.5, length(x)))
  (x - lo) / (hi - lo)
}

dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(args[[2]], width = 9, height = 9, family = "sans")
circlize::circos.clear()
circlize::circos.par(start.degree = 90, gap.after = rep(2, length(chrom)),
                     cell.padding = c(0, 0, 0, 0), track.margin = c(0.006, 0.006))
circlize::circos.initialize(factors = chrom,
                            xlim = cbind(rep(0, length(chrom)), unname(lengths)))
# Karyotype ring: one sector per chromosome, with labels in the same order as
# the supplied chromosome-size table / window file.
circlize::circos.trackPlotRegion(factors = chrom, ylim = c(0, 1),
  track.height = 0.055, bg.col = "#DADDE0", bg.border = "white",
  panel.fun = function(x, y) {
    sector <- circlize::CELL_META$sector.index
    circlize::circos.text(lengths[[sector]] / 2, 0.5, sector,
                          facing = "bending.inside", cex = 0.42, niceFacing = TRUE)
  })
for (metric in tracks) {
  values <- rescale(as.numeric(d[[metric]]), metric)
  d$scaled <- values
  circlize::circos.trackPlotRegion(factors = chrom, ylim = c(0, 1),
    track.height = 0.11, bg.border = "#E0E0E0",
    panel.fun = function(x, y) {
      sector <- circlize::CELL_META$sector.index
      z <- d[d$chrom == sector, , drop = FALSE]
      xx <- (z$start + z$end) / 2
      ok <- is.finite(z$scaled)
      if (any(ok)) circlize::circos.lines(xx[ok], z$scaled[ok],
                                          col = colors[[metric]], lwd = 0.75)
    })
}
graphics::title(main = args[[3]], cex.main = 1.1)
graphics::legend("center", legend = c("GERP", "Genes", "TE", "SNP", "Nucleotide diversity"),
                 col = colors, lty = 1, lwd = 2, bty = "n", cex = 0.75)
circlize::circos.clear()
grDevices::dev.off()
