# Heterozygosity, ROH and burden

`selected/common/calculate_pi_heterozygosity.sh` calculates π and heterozygosity. The species-specific `run_roh_froh.sh` scripts calculate ROH and FROH with PLINK gap100/kb100 and five missing sites per window.

`plot_heterozygosity_roh_fig6_abc.R` plots Fig. 6A–C. `plot_bamboo_roh_landscape_fig6d.R` plots the bamboo-only Fig. 6D. For Supplementary Fig. S8, `prepare_s8_source.py` checks the burden data against Supplementary Table S6 and `plot_s8_genotype_burden.R` draws the panels. Supply the final S3/S6 workbooks and matching analysis tables locally.
