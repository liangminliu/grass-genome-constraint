# Whole-genome alignment

`selected/run_msa_with_existing_tree.sh` runs the reference-guided alignment using an existing guide tree. The four original `call_conservation.smk` rules and their shared helper scripts are preserved in `selected/msa_pipeline_source/`; they contain the trimAl, RAxML, phyloFit and gerpcol commands used by the Snakemake workflow. Use the matching focal-reference rule within the original msa_pipeline checkout. `selected/coverage/` contains Fig. 2 coverage and alignment-depth calculations.
