# Figure 2A: observed MAF coverage by query species

This subflow runs once per reference-coordinate alignment (`ped`, `Zmay`, or `Atau`). It produces one coverage plot per reference; assembling the three panels into the manuscript Fig. 2A is a separate layout step. Fig. 2B's cumulative alignment-depth calculation is separate.

1. `maf_pair_coverage.sbatch` runs `mafPairCoverage --bin_length 5000` for one query species per array task. Supply `REF_SPECIES`, `MAF_FILE`, `SPECIES_LIST` and optionally `OUT_DIR` through the job environment. Submit with `--array=1-N`, where N is the number of query rows; create `logs/` before submission.
2. `summarize_coverage.py` reads the reference row (`ped*`, `Zmay*`, or `Atau*`) from the sixth field of each result, converts the fractional value to percent and writes `Species`/`ObsCoverage` in the supplied query order. It requires one valid result per query.
3. `plot_grouped_coverage_final.R` joins that table to the two-column phylogenetic-order/subfamily list and writes a PDF. The original eight-color palette and compact bar-chart layout are retained. The y-axis expands beyond 80% when needed so bars are not dropped.

Example for bamboo, using the supplied query classification list:

```bash
CODE_ROOT=/path/to/grass-genome-constraint
CLASS_LIST="$CODE_ROOT/metadata/94sp.pop.list"
TOOLS="$CODE_ROOT/scripts/02_alignment/selected/coverage"
mkdir -p logs mafPairCoverage_results
N=$(wc -l < "$CLASS_LIST")  # 93 queries, excluding the bamboo reference
sbatch --array=1-${N} --export=ALL,REF_SPECIES=ped,MAF_FILE=roast.maf,SPECIES_LIST="$CLASS_LIST" "$TOOLS/maf_pair_coverage.sbatch"
python3 "$TOOLS/summarize_coverage.py" --reference ped --species-list "$CLASS_LIST" --results-dir mafPairCoverage_results --output sorted_coverage.tsv
Rscript "$TOOLS/plot_grouped_coverage_final.R" "$CLASS_LIST" sorted_coverage.tsv final_grouped_coverage.pdf
```

For maize and goatgrass, use the same scripts with their own MAF, reference ID and query list. Each query list must exclude its own reference and include the appropriate bamboo representation; the supplied bamboo list cannot simply be reused because it contains `Zmay` and `Atau` but no bamboo reference. Confirm the MAF sequence prefixes and final 94-genome panel before constructing those lists.

The supplied `94sp.pop.list` has **93 distinct rows**. In the bamboo analysis, adding the reference gives 94 grass genomes; the historical `--array=1-94` requested one extra task. Historical versions are withheld from this confirmed release.

## Figure 2B: cumulative alignment depth

`compute_alignment_depth.py` reads chromosome BED files with one header row and aligned-species count `k` in column 6. It writes a `k`/`All` cumulative percentage table for one reference; run the same program for bamboo, maize and goatgrass. By default each BED interval contributes its length in bases, consistent with the manuscript's proportion of eligible reference sites. `--weight rows` reproduces the author-supplied `all.py` calculation, which weights every BED record equally. The two calculations agree when every BED interval is one base long. The program reports how many intervals are not one base; record this before choosing the final figure source.

`plot_alignment_depth_three_species.py` joins the three resulting curves into one Fig. 2B PDF. Earlier fixed-path and row-weighted scripts are withheld from this confirmed release.

```bash
python3 "$TOOLS/compute_alignment_depth.py" --bed-dir PED_BED_DIR --output ped_depth.tsv --max-k 94
python3 "$TOOLS/compute_alignment_depth.py" --bed-dir ZMAY_BED_DIR --output zmay_depth.tsv --max-k 94
python3 "$TOOLS/compute_alignment_depth.py" --bed-dir ATAU_BED_DIR --output atau_depth.tsv --max-k 94
python3 "$TOOLS/plot_alignment_depth_three_species.py" --ped ped_depth.tsv --zmay zmay_depth.tsv --atau atau_depth.tsv --output figure2b_alignment_depth.pdf
```

The BEDs must represent the same eligible reference-site definition as the Fig. 2B denominator. Their source folders and final hashes have not been confirmed. If the available BEDs contain only GERP-scored sites or overlapping intervals, the resulting curve may have a different denominator; check the inputs against the final figure before treating this as its exact generating run.
