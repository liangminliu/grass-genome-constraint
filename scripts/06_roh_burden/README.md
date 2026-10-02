# π, heterozygosity, ROH and FROH

The three `selected/{ped,teosinte,atau}/pi_het.sh` jobs use `selected/common/calculate_pi_heterozygosity.sh` to calculate site π, mean π and individual heterozygosity from each species' filtered VCF. Create `logs/` before SLURM submission.

The three `run_roh_froh.sh` scripts are the author-supplied final ROH/FROH calculations. They execute PLINK gap100/kb100 with `--homozyg-window-missing 5`. Supply the VCF/PLINK inputs and each species' `Lauto.from_filtered_SNPs.tsv` as described in those scripts. A. tauschii's ROH job excludes `BW_01192` from the 141-ID input.

## Heterozygosity and ROH visualization (Results section 5; Fig. 6A-C)

`plot_heterozygosity_roh_fig6_abc.R` recreates the current manuscript's revised Fig. 6A-C directly from Supplementary Table S3B-D. It plots individual observed heterozygosity, FROH distributions, and cumulative FROH by four ROH length classes for the final 63 bamboo, 90 teosinte, and 140 A. tauschii individuals. It checks sample counts, duplicate IDs, numeric ranges, and whether the four class fractions sum to each individual's FROH. It does not remove an extreme bamboo sample or truncate the FROH axis.

```bash
Rscript scripts/06_roh_burden/plot_heterozygosity_roh_fig6_abc.R \
  /path/to/Supplementary_Table_S3_Population_Samples_Heterozygosity_ROH.xlsx \
  results/Fig6_heterozygosity_ROH_ABC
```

The output is a standalone PDF and PNG. R `readxl` is required. The source workbook is not included because it contains individual-level study data; obtain the approved S3 workbook before reproduction. The current Fig. 6D chromosome track was retained from an earlier composition and lacks a confirmed standalone final plotting command and representative-sample input. Those are still requested. The older species-specific ROH plotting candidates are therefore not presented as the final Fig. 6 source.

Individual GERP/SIFT burden and final Fig. S8 table assembly are pending exact production files and are not included here.
