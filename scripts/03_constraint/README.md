# Evolutionary constraint

`selected/region_annotation/` intersects GERP thresholds with genomic categories and plots regional proportions. `selected/gerp_histogram/` plots the bamboo GERP distribution for Fig. 4A.

For Fig. 3/S1, run `selected/window_tracks/compute_500kb_tracks.py` and `plot_circular_tracks.R` for each focal genome. These 500-kb window scripts implement the Methods definitions and require matching GERP, gene, TE, SNP and nucleotide-diversity inputs. BED and bedGraph coordinates must be zero-based, half-open.
