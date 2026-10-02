# Result-figure code audit

The code below was checked against the current manuscript figure captions and the available script inventories. A historical figure filename alone does not establish that its code produced the current panel.

| Current panel | Released code | Remaining source or check |
| --- | --- | --- |
| Fig. 1B | `scripts/01_synteny/plot_genetribe_three_references.R` | Verify the final phylogeny order and 13-genome chromosome view (Fig. 1C). |
| Fig. 2A-B | `scripts/02_alignment/` coverage and alignment-depth plots | Supply the final three MAFs and panel input tables. |
| Fig. 3; S1 | `scripts/03_constraint/selected/window_tracks/` now computes and plots all 500-kb tracks from Methods | This is a Methods reconstruction; supply original per-base/annotation/VCF inputs and compare its output to the current image. |
| Fig. 4 | `scripts/03_constraint/` contains supplied GERP distribution and compartment plots | Confirm exact current panel-to-output correspondence for B-D. |
| Fig. 5 | Historical allele-frequency, annotation and bamboo burden plots found | Current four-panel assembly and input-to-panel mapping are not yet verified; these candidates were not promoted. |
| Fig. 6A-C | `scripts/06_roh_burden/plot_heterozygosity_roh_fig6_abc.R` | Approved S3B-D workbook required. |
| Fig. 6D | `scripts/06_roh_burden/plot_bamboo_roh_landscape_fig6d.R` | Bamboo FROH, chromosome and ROH segment tables required; this panel is bamboo only. |
| Fig. 7; S9; S10 | `scripts/07_gene_compare/selected/` now includes the retained strict-core and gene-class analysis/plotting programs plus panel-to-test audit | 5,108 groups, 668 matched SNP-positive groups and GERP >2 maize–A. tauschii adjusted P=0.097203 were verified; inspect the existing image brackets and black markers. |
| S8 | `scripts/06_roh_burden/prepare_s8_source.py` and `plot_s8_genotype_burden.R` | Reproduces 293 final source rows from approved S6/raw tables; rerun the R plot in an environment with R. |

For Fig. 6D, the executed bamboo script named the corresponding representative landscape its local panel E. The current manuscript labels that bamboo-only track Fig. 6D. The released code retains the original five-high/five-low non-`PE` display rule while validating membership against the final 63-individual bamboo set.
