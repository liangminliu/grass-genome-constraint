# Constrained-site region summaries and plots

1. `run_intersect.sh`: SLURM array tasks 1–8 cover PedC/PedD and GERP >0, >2, >4, >6 BED files. Each task intersects one header-bearing GERP BED with the category BEDs and writes per-category overlap files and `category_stats.tsv`.
2. Prepare the curated `gerp_constrained_master.csv` table from the region statistics and final category definitions. Its production script remains to be supplied.
3. `plot_CD_combined.R`: reads `gerp_constrained_master.csv` and writes `GERP_pie_C_D_combined.pdf` for the bamboo C/D implementation detail.
4. `plot_gt2_multi_species.R`: reads a headerless three-column `species`, `region`, `count` TSV (default `all_gt2_CDS.tsv`) and writes `GERP_gt2_species_pies.pdf`. Run `Rscript plot_gt2_multi_species.R [input.tsv] [output.pdf]`. It checks species/region identifiers, duplicate pairs and non-negative counts, fills absent regions with zero, and collects one shared legend. The script that assembled the exact three-species input has not been identified.

The `-wa -wb` intersection counts **overlap rows**, not unique GERP bases. A GERP interval overlapping multiple category intervals can contribute multiple rows, and categories may overlap one another. `proportion` is each category's fraction of the sum of these overlap rows. Before treating the pies as mutually exclusive genomic-region fractions, verify disjoint category BEDs or establish an explicit category-priority assignment and regenerate the counts. `GERP>0` is retained as a source-script exploratory threshold; the manuscript's cumulative thresholds are >2, >4 and >6.

These plots are code provenance for the bamboo C/D and three-species summaries. They do not by themselves establish that the PDFs are the exact manuscript panels. Bamboo is interpreted as one whole-genome focal species in the main Results.

The three-species plot expects the source labels `P.edulis`, `Z.may`, and `A.tauschii`; the maize label is kept exactly as in the supplied plotting code. Its five region labels are `CDS`, `Intron`, `UTR`, `UpDown5K`, and `Inter-genic`. Confirm that the input table uses these labels and an appropriate, consistent site-counting rule before drawing the manuscript figure.
