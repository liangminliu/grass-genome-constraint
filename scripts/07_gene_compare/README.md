# Synteny-based gene comparisons

The confirmed SynPan run uses five labels, in order: `Atau`, `osat`, `pedC`, `pedD`, `Zmay`. Run `selected/synpan/run_synpan.sbatch` with matching `<label>.prot` and `<label>.bed` inputs; see `selected/synpan/README.md`. The main `selected/synpan_build.pl` and 14 companion Perl scripts are included.

## Fig. 7 and Supplementary Figs. S9-S10

The source inventory now includes two full analysis/plotting programs under `selected/`:

- `plot_strict_core_fig7_s9.R` pools the two bamboo homeologs per strict-core group, builds gene/group tables, tests matched groups and draws strict-core comparisons. The source script has an internal `v10` identifier and writes both tables and figures.
- `plot_gene_classes_fig7_s10.R` builds the within-species strict-core/relaxed-syntenic/nonsyntenic gene-class tables, Kruskal-Wallis and BH-adjusted pairwise Wilcoxon results, and figures.

Both scripts expect the source files beside the working directory, including `core_5sp_1to1.SG.pan`, GERP CDS-per-gene tables, SNP/deleterious counts and the gene-class lists. Run them from a copy of that input directory; their default `OUTDIR` values are local relative paths. Large gene-level inputs are not committed.

`selected/audit_fig7_sources.py` validates the resulting tables and writes `panel_to_source_test.tsv` and `source_audit.json`. It was run against the retained source tables: 5,108 complete strict-core groups at GERP >4, 668 matched SNP-positive groups, the two documented A. tauschii zero-filled SNP records (SG0039462 and SG0039466), and 12 required table/statistic files were confirmed. The GERP >2 Z. mays–A. tauschii BH-adjusted P value is **0.0972031305**, so that S9 comparison must not be marked significant. The checked mapping and audit summary are committed under `selected/`.

```bash
python scripts/07_gene_compare/selected/audit_fig7_sources.py \
  --strict-dir /path/to/strict_core_results \
  --class-dir /path/to/gene_class_results \
  --outdir results/fig7_audit
```

The plotted significance brackets and the meaning of black markers in the existing manuscript image still require a source-image check. The source scripts and test tables are now recoverable, but the existing assembled Fig. 7/S9/S10 images should not be described as fully validated until that visual check is complete.
