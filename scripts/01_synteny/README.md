# Genome-wide synteny: GeneTribe

This directory contains one consolidated GeneTribe launcher for three focal references, the exact 94-name order supplied by the author, a result summarizer, and the renamed grouped plotting script. The launcher replaces the historical bamboo-only 81-species batch and later two-reference chain. It schedules **94 query jobs, 279 non-self comparisons, and 93 comparisons per reference**.

| File | Role |
| --- | --- |
| `species_94.list` | Exact author-supplied names and order; case is preserved. |
| `run_genetribe_three_references.sh` | Validate inputs, generate a pair manifest, or run one query against all non-self references. |
| `summarize_genetribe_three_references.py` | Summarize directional blocks, RBH, one-to-many and singleton outputs; self-comparisons are `NA`. |
| `plot_genetribe_three_references.R` | Produce grouped three-reference syntenic-gene and gene-composition plots. |

Place `<species>.fa`, `<species>.bed`, and `<species>.chrlist` for every exact name in `species_94.list` in one data directory. Provide a separate two-column `genome.size` table (`species`, genome size in bp) for the summary. Install GeneTribe in the runtime environment. The genome data and executable are not included.

```bash
export DATA_DIR=/path/to/data
export WORK_DIR=/path/to/genetribe_results
bash scripts/01_synteny/run_genetribe_three_references.sh validate
bash scripts/01_synteny/run_genetribe_three_references.sh manifest > genetribe_pairs.tsv
mkdir -p logs
bsub -J 'genetribe[1-94]' -q Q96C1T_X12 -n 8 \
  -o 'logs/genetribe.%J.%I.out' -e 'logs/genetribe.%J.%I.err' \
  'bash scripts/01_synteny/run_genetribe_three_references.sh run'
```

The launcher stores each pair in `<query>/ref_<reference>_output/`, skips completed pairs, and stops if an incomplete result already exists. It uses the observed orientation `genetribe core -l <reference> -f <query>`. Three reference runs are sequential within each query job to avoid GeneTribe temporary-file collisions. Adapt the LSF queue and core request to the cluster.

After all array jobs finish:

```bash
cd "$WORK_DIR"
python /path/to/repo/scripts/01_synteny/summarize_genetribe_three_references.py \
  -l /path/to/repo/scripts/01_synteny/species_94.list \
  -g /path/to/genome.size -d "$DATA_DIR" -o genetribe_main_summary.tsv
Rscript /path/to/repo/scripts/01_synteny/plot_genetribe_three_references.R \
  genetribe_main_summary.tsv \
  /path/to/repo/scripts/01_synteny/species_94.list \
  genetribe_three_reference_plots
```

The R script requires `ggplot2`, `dplyr`, `tidyr`, `scales`, and `patchwork`. It retains empty species slots and omits self-comparison bars. Inspect the summary `status` column before using the plots: only `OK` pairs contribute values. The final collinear-block retention rule and JCVI chromosome-view command are still awaiting confirmation.
