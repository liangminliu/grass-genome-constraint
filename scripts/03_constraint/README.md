# GERP-based evolutionary constraint

`selected/region_annotation/` overlaps bamboo C/D cumulative-threshold BEDs with genomic regions and plots bamboo and three-species compositions. `selected/gerp_histogram/` plots the author-supplied bamboo GERP score distribution used for Fig. 4A. The scoring pipeline and neutral model inputs are held until their exact production versions are supplied; the released scripts begin with GERP output files.

Classes are cumulative GERP >2, >4 and >6. Region pies require category rules that make fractions interpretable; supply the final category BEDs and source tables.

## Fig. 3 and Fig. S1 genome tracks

`selected/window_tracks/compute_500kb_tracks.py` builds the non-overlapping 500-kb windows described in Methods for bamboo (Fig. 3A), maize (Fig. 3B), and A. tauschii (Fig. S1). `selected/window_tracks/plot_circular_tracks.R` draws the chromosome, GERP, gene, TE, SNP, and nucleotide-diversity rings from each table. Use the same two scripts for all three species.

```bash
python scripts/03_constraint/selected/window_tracks/compute_500kb_tracks.py \
  --chrom-sizes species.chrom.sizes --gerp-bedgraph species.gerp.bedgraph.gz \
  --genes-bed species.genes.bed --te-bed species.te.bed \
  --snp-vcf species.filtered.snps.vcf.gz --pi-bedgraph species.pi.bedgraph.gz \
  --output results/species.500kb_tracks.tsv
Rscript scripts/03_constraint/selected/window_tracks/plot_circular_tracks.R \
  results/species.500kb_tracks.tsv results/species.genome_tracks.pdf "Species name"
```

All BED/bedGraph intervals must be zero-based half-open and use the same chromosome names and assembly as the VCF (whose POS field is one-based). GERP and π tracks are means weighted by scored base pairs; a window with no scored bases is `NA`, not zero. Gene and TE tracks are the fraction of window bases covered by the merged respective annotations. SNP density is unique biallelic SNP records per kb of physical window length. Input GERP/π bedGraphs must be sorted and non-overlapping per chromosome. The R plot requires `circlize`. These scripts implement the Methods definitions; the exact source annotations and plotted outputs must be checked against the existing manuscript image before labeling them a pixel-identical reconstruction.
