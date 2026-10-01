#!/usr/bin/env bash
#SBATCH --job-name=atau_functional_intersect
#SBATCH --cpus-per-task=8
#SBATCH --array=1-3%3
#SBATCH --output=logs/atau_functional_%A_%a.out
#SBATCH --error=logs/atau_functional_%A_%a.err

# Summarize GERP-positive SNP positions across functional-region BEDs.
# Array size is three because the inputs are merged_gt2, gt4 and gt6.
set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit as an array job}"
thresholds=(merged_gt2 merged_gt4 merged_gt6)
idx=$((SLURM_ARRAY_TASK_ID - 1))
[[ $idx -ge 0 && $idx -lt ${#thresholds[@]} ]] || exit 2
name=${thresholds[$idx]}
input="$name/${name}_gerp_regions.bed"
feature_dir=${FEATURE_DIR:-genomic_features_results}
[[ -s "$input" && -d "$feature_dir" ]] || {
  echo "Missing GERP regions or functional BED directory" >&2; exit 1;
}
command -v bedtools >/dev/null || { echo "bedtools is required" >&2; exit 1; }
shopt -s nullglob
features=("$feature_dir"/*processed.bed)
(( ${#features[@]} > 0 )) || { echo "No processed functional BEDs" >&2; exit 1; }

outdir="$name/functional_analysis"
mkdir -p "$outdir" logs
summary="$outdir/${name}_summary.tsv"
printf 'Region\tOverlap_intervals\tGERP_regions\tFUNC_regions\tOverlap_bp\tGERP_coverage(%%)\tFUNC_coverage(%%)\n' > "$summary"

# Use union lengths for coverage denominators; retain detailed overlap rows
# separately so the original interval-level outputs remain inspectable.
gerp_total=$(bedtools sort -i "$input" | bedtools merge -i stdin |
  awk '{sum += $3 - $2} END {print sum+0}')
for feature in "${features[@]}"; do
  region=$(basename "$feature" _processed.bed)
  overlap="$outdir/${name}_${region}_overlap.bed"
  gerp_out="$outdir/${name}_${region}_gerp_regions.bed"
  feature_out="$outdir/${name}_${region}_func_regions.bed"
  detail="$outdir/${name}_${region}_detailed.bed"
  bedtools intersect -a "$input" -b "$feature" -wa -u > "$gerp_out"
  bedtools intersect -a "$feature" -b "$input" -wa -u > "$feature_out"
  bedtools intersect -a "$input" -b "$feature" -wa -wb > "$detail"
  bedtools intersect -a "$input" -b "$feature" |
    bedtools sort -i stdin | bedtools merge -i stdin > "$overlap"
  overlap_n=$(wc -l < "$overlap")
  gerp_n=$(wc -l < "$gerp_out")
  feature_n=$(wc -l < "$feature_out")
  overlap_bp=$(awk '{sum += $3 - $2} END {print sum+0}' "$overlap")
  feature_total=$(bedtools sort -i "$feature" | bedtools merge -i stdin |
    awk '{sum += $3 - $2} END {print sum+0}')
  gerp_pct=$(awk -v n="$overlap_bp" -v d="$gerp_total" 'BEGIN {printf "%.4f", d ? 100*n/d : 0}')
  feature_pct=$(awk -v n="$overlap_bp" -v d="$feature_total" 'BEGIN {printf "%.4f", d ? 100*n/d : 0}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$region" "$overlap_n" "$gerp_n" "$feature_n" "$overlap_bp" "$gerp_pct" "$feature_pct" >> "$summary"
done
echo "Wrote $summary"
