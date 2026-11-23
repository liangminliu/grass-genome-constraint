# Fig. 2 coverage

`maf_pair_coverage.sbatch` computes pairwise MAF coverage; `summarize_coverage.py` and `plot_grouped_coverage_final.R` produce the grouped coverage plot. `compute_alignment_depth.py` and `plot_alignment_depth_three_species.py` produce cumulative alignment-depth curves. Run the workflow separately for each focal reference with its matching MAF and species list.

For Supplementary Fig. S2, run plot_s2_coverage_distance.R from the directory containing the three reference-specific MAF coverage folders and neutral trees. ObsCoverage_tripanel.pdf is the manuscript panel; the script also writes a Coverage_tripanel.pdf diagnostic and source/statistics TSV files.
