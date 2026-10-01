# Whole-genome alignment and Fig. 2 coverage

`selected/run_msa_with_existing_tree.sh` uses a supplied ROADIES tree with the original MSA pipeline checkout. `selected/prepare_existing_tree_config.py` prunes that guide topology to the configured species. The runner orders `align`, `roast` and `call_conservation`; it requires the original `workflow/Snakefile` and focal config YAML, which are not bundled until their final revisions are identified.

`selected/coverage/` holds the author-supplied coverage job, summary and plots for Fig. 2A, plus cumulative alignment-depth calculations/plots for Fig. 2B. The 93 query-species metadata rows are in `metadata/94sp.pop.list`. Supply the final three reference MAF/BED inputs and check whether BED rows represent bases or segments before interpreting the depth denominator.
