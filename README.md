# Evolutionary constraint and deleterious variation across grass genomes

Scripts associated with Liu et al., *Deep evolutionary constraint predicts deleterious variation in grass genomes*. The workflow follows the Methods in `MANUSCRIPT_MASTER_20260929_abstract_references_tracked.docx`.

| Step | Code |
| --- | --- |
| Grass phylogeny | [`scripts/00_phylogeny/`](scripts/00_phylogeny/) |
| Genome-wide synteny | [`scripts/01_synteny/`](scripts/01_synteny/) |
| Whole-genome alignment | [`scripts/02_alignment/`](scripts/02_alignment/) |
| Evolutionary constraint | [`scripts/03_constraint/`](scripts/03_constraint/) |
| Population SNP filtering | [`scripts/04_population_snps/`](scripts/04_population_snps/) |
| Deleterious variant annotation | [`scripts/05_annotation/`](scripts/05_annotation/) |
| ROH, FROH and individual burden | [`scripts/06_roh_burden/`](scripts/06_roh_burden/) |
| Synteny-based gene comparisons | [`scripts/07_gene_compare/`](scripts/07_gene_compare/) |

The directories contain analysis and plotting scripts for Figs. 1–7 and their supplementary figures. Input genomes, sequence alignments, VCFs and individual-level results are not included. Set the paths in the scripts for the target computing environment.

Key settings: ROADIES continuity 85 and support threshold 0.95; GeneTribe compares each focal reference with 93 other grasses; GERP thresholds >2, >4 and >6; SIFT <0.05; PLINK ROH gap100/kb100 with five missing sites per window.

Some bamboo and maize upstream SNP commands and the Fig. 3/S1 window-track scripts were reconstructed from Methods. Their directory READMEs identify them. The remaining scripts were consolidated from supplied analysis code.

## Citation

Liu L-M et al. *Deep evolutionary constraint predicts deleterious variation in grass genomes*. Manuscript in preparation.
