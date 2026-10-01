#!/usr/bin/env bash
#SBATCH --job-name=gerp_intersect
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-8
#SBATCH --output=logs/intersect_%A_%a.out
#SBATCH --error=logs/intersect_%A_%a.err

# Cleaned release version of run_intersect_original.sh. The historical count
# is the number of -wa -wb overlap rows, not unique GERP sites.
set -euo pipefail
module load bedtools 2>/dev/null || true
command -v bedtools >/dev/null || { echo "bedtools is required" >&2; exit 1; }

GERP_DIR=${GERP_DIR:-/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/84sp_ped/94sp_ped_results/GERP_final_output}
CAT_DIR=${CAT_DIR:-../2.1.constrained_category}
OUT_BASE=${OUT_BASE:-./intersect_results}
task=${SLURM_ARRAY_TASK_ID:-${1:-}}
[[ "$task" =~ ^[1-8]$ ]] || { echo "Specify array task 1-8" >&2; exit 2; }

files=(
  "$GERP_DIR/gerp_gt0/pedC_gt0.bed" "$GERP_DIR/gerp_gt0/pedD_gt0.bed"
  "$GERP_DIR/gerp_gt2/pedC_gt2.bed" "$GERP_DIR/gerp_gt2/pedD_gt2.bed"
  "$GERP_DIR/gerp_gt4/pedC_gt4.bed" "$GERP_DIR/gerp_gt4/pedD_gt4.bed"
  "$GERP_DIR/gerp_gt6/pedC_gt6.bed" "$GERP_DIR/gerp_gt6/pedD_gt6.bed"
)
gerp=${files[task-1]}
[[ -s "$gerp" ]] || { echo "Missing GERP BED: $gerp" >&2; exit 1; }
shopt -s nullglob
categories=("$CAT_DIR"/*.bed)
(( ${#categories[@]} > 0 )) || { echo "No category BEDs in $CAT_DIR" >&2; exit 1; }

group=$(basename "$(dirname "$gerp")")
base=$(basename "$gerp" .bed)
out_dir="$OUT_BASE/$group"
mkdir -p "$out_dir"
clean=$(mktemp "$out_dir/${base}.clean.XXXXXX")
trap 'rm -f "$clean"' EXIT
tail -n +2 "$gerp" > "$clean"  # source GERP BED has one header row
[[ -s "$clean" ]] || { echo "GERP BED has no data rows: $gerp" >&2; exit 1; }

stats="$out_dir/${base}.category_stats.tsv"
counts="$out_dir/${base}.counts.tmp"
: > "$counts"
pids=()
for category in "${categories[@]}"; do
  name=$(basename "$category" .bed)
  [[ -s "$category" ]] || { echo "Empty category BED: $category" >&2; exit 1; }
  output="$out_dir/${base}.${name}.intersect.bed"
  (
    bedtools intersect -wa -wb -a "$clean" -b "$category" > "$output"
    printf '%s\t%s\n' "$name" "$(wc -l < "$output")" > "$out_dir/${base}.${name}.count.tmp"
  ) &
  pids+=("$!")
  if (( ${#pids[@]} >= ${SLURM_CPUS_PER_TASK:-10} )); then
    for pid in "${pids[@]}"; do wait "$pid"; done
    pids=()
  fi
done
for pid in "${pids[@]}"; do wait "$pid"; done

for category in "${categories[@]}"; do
  name=$(basename "$category" .bed)
  cat "$out_dir/${base}.${name}.count.tmp" >> "$counts"
  rm -f "$out_dir/${base}.${name}.count.tmp"
done
total=$(awk -F '\t' '{s+=$2} END {print s+0}' "$counts")
{
  printf 'category\tcount\tproportion\n'
  awk -F '\t' -v total="$total" 'BEGIN{OFS="\t"} {p=(total ? $2/total : 0); printf "%s\t%d\t%.6f\n",$1,$2,p}' "$counts"
} > "$stats"
rm -f "$counts"
echo "Wrote $stats; overlap rows across categories: $total"
