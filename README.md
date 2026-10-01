# Evolutionary constraint and deleterious variation across grass genomes

Code associated with Liu et al., *Deep evolutionary constraint predicts deleterious variation in grass genomes*.

This repository follows the Methods order. It contains the author-supplied scripts and the compact workflows designated for release. Historical alternatives, exploratory reruns and methods-only drafts are held outside this repository. `CONFIRMED_FILES.tsv` records the source and SHA-256 of every included code or metadata file. Large genome, MAF and VCF inputs are not included.

| Step | Directory | Released workflow |
| --- | --- | --- |
| Grass phylogeny | [`scripts/00_phylogeny/`](scripts/00_phylogeny/) | ROADIES accurate mode, 97 input genomes |
| 1. Genome-wide synteny | [`scripts/01_synteny/`](scripts/01_synteny/) | 3 focal references × 93 GeneTribe queries; final job manifest requested |
| 2. Whole-genome alignment | [`scripts/02_alignment/`](scripts/02_alignment/) | Existing-tree MSA driver and Fig. 2 coverage/depth scripts |
| 3. GERP-based evolutionary constraint | [`scripts/03_constraint/`](scripts/03_constraint/) | GERP category intersections and score figures; scoring inputs requested |
| 4. Population SNP filtering | [`scripts/04_snp_filter/`](scripts/04_snp_filter/) | A. tauschii filtered VCF and exact 141-name retained list |
| 5. Deleterious variant annotation | [`scripts/05_annotation/`](scripts/05_annotation/) | Three-species SIFT 4G and SnpEff workflows; A. tauschii SNP–GERP overlap |
| 6. ROH, FROH and individual burden | [`scripts/06_roh_burden/`](scripts/06_roh_burden/) | Three-species π/heterozygosity and author-supplied final ROH/FROH scripts |
| 7. Synteny-based gene comparisons | [`scripts/07_gene_compare/`](scripts/07_gene_compare/) | Five-input SynPan and its 14 companion Perl scripts |

Read each step's `README.md` before running it. The steps currently missing a confirmed production script or exact input are listed in [`MISSING_INPUTS_AND_CODE.md`](MISSING_INPUTS_AND_CODE.md); a Methods mention alone is not treated as evidence of an executed command.

## Confirmed analysis settings

- ROADIES continuity 85, support threshold 0.95 and initial `GENE_COUNT: 16000` sampled fragments.
- GeneTribe: each focal reference compared with 93 query grass genomes.
- GERP classes: cumulative >2, >4 and >6. SIFT damaging prediction: <0.05.
- PLINK ROH: `--homozyg-gap 100`, `--homozyg-kb 100`, `--homozyg-window-missing 5`.
- The supplied A. tauschii retention list has 141 unique IDs, including `BW_01192`; the ROH script excludes that ID for its 140-individual analysis.
- SynPan labels: `Atau`, `osat`, `pedC`, `pedD`, `Zmay`. Bamboo C/D are two subgenome inputs, not separate focal species.

Scripts retain HPC paths from the source runs. Set cluster paths and input names as described in the step README files before submission. SIFT/GERP labels are predictions, not measured fitness effects.

## Manuscript figures

| Figure | Analysis |
| --- | --- |
| Fig. 1 | Phylogeny and synteny |
| Fig. 2 | Whole-genome alignment and coverage |
| Figs. 3–4 | Evolutionary constraint; the supplied GERP distribution is Fig. 4A |
| Fig. 5 | Predicted deleterious variation |
| Fig. 6 | Heterozygosity and ROH |
| Fig. 7 | Synteny-based constraint and coding variation |

Major software includes ROADIES, treePL, GeneTribe, JCVI, LAST, MULTIZ/ROAST, PHAST, GERP++, BWA-MEM, GATK, VCFtools, PLINK, SIFT 4G and SnpEff. The repository contains confirmed code for only the steps identified in the table above; [`MISSING_INPUTS_AND_CODE.md`](MISSING_INPUTS_AND_CODE.md) identifies the remaining production commands and data.

## Citation

Liu L-M et al. *Deep evolutionary constraint predicts deleterious variation in grass genomes*. Manuscript in preparation.
