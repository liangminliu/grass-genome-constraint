# π, heterozygosity, ROH and FROH

The three `selected/{ped,teosinte,atau}/pi_het.sh` jobs use `selected/common/calculate_pi_heterozygosity.sh` to calculate site π, mean π and individual heterozygosity from each species' filtered VCF. Create `logs/` before SLURM submission.

The three `run_roh_froh.sh` scripts are the author-supplied final ROH/FROH calculations. They execute PLINK gap100/kb100 with `--homozyg-window-missing 5`. Supply the VCF/PLINK inputs and each species' `Lauto.from_filtered_SNPs.tsv` as described in those scripts. A. tauschii's ROH job excludes `BW_01192` from the 141-ID input.

Individual GERP/SIFT burden and final Fig. 6/S8 table assembly are pending exact production files and are not included here.
