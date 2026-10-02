# Genome-wide synteny

`run_genetribe_three_references.sh` runs GeneTribe for the three focal references against the other 93 grasses in `species_94.list`. `summarize_genetribe_three_references.py` combines the results, and `plot_genetribe_three_references.R` makes the grouped plots. Supply the genome FASTA, gene BED and chromosome list for each species; set `DATA_DIR` and `WORK_DIR` before running the LSF array.
