#!/usr/bin/env bash
# One production entry point for the three GeneTribe reference comparisons.
# Submit one array task per query: bsub -J 'genetribe[1-94]' -n 8 -q Q96C1T_X12 \
#   -o 'logs/genetribe.%J.%I.out' -e 'logs/genetribe.%J.%I.err' \
#   'bash scripts/01_synteny/run_genetribe_three_references.sh run'
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPECIES_LIST="${SPECIES_LIST:-${SCRIPT_DIR}/species_94.list}"
DATA_DIR="${DATA_DIR:-${PWD}/data}"
WORK_DIR="${WORK_DIR:-${PWD}/genetribe_results}"
REFS=(Phyllostachys_edulis Aegilops_tauschii Zea_mays)
SUFFIXES=(block_pos collinearity_info one2many one2one SBH RBH singleton)

usage() {
  echo "Usage: DATA_DIR=/path/to/data WORK_DIR=/path/to/results $0 validate|manifest|run [query_index]" >&2
  exit 2
}

[[ $# -ge 1 ]] || usage
mode="$1"
mapfile -t species < <(sed '/^[[:space:]]*$/d' "$SPECIES_LIST")
[[ ${#species[@]} -eq 94 ]] || { echo "Expected 94 species, found ${#species[@]}" >&2; exit 1; }
declare -A seen=()
for sp in "${species[@]}"; do
  [[ "$sp" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Invalid species name: $sp" >&2; exit 1; }
  key="${sp,,}"
  [[ ! ${seen[$key]+present} ]] || { echo "Duplicate species name: $sp" >&2; exit 1; }
  seen[$key]=1
done

for ref in "${REFS[@]}"; do
  [[ ${seen[${ref,,}]+present} ]] || { echo "Reference absent from species list: $ref" >&2; exit 1; }
done

if [[ "$mode" == manifest ]]; then
  printf 'query_index\tquery_species\treference_species\n'
  for i in "${!species[@]}"; do
    for ref in "${REFS[@]}"; do
      [[ "${species[$i],,}" == "${ref,,}" ]] || printf '%d\t%s\t%s\n' "$((i+1))" "${species[$i]}" "$ref"
    done
  done
  exit 0
fi

[[ "$mode" == validate || "$mode" == run ]] || usage
[[ -d "$DATA_DIR" ]] || { echo "Missing data directory: $DATA_DIR" >&2; exit 1; }
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
for sp in "${species[@]}"; do
  for ext in fa bed chrlist; do
    [[ -s "$DATA_DIR/$sp.$ext" ]] || { echo "Missing or empty: $DATA_DIR/$sp.$ext" >&2; exit 1; }
  done
done
echo "Validated 94 species and 282 input files."
[[ "$mode" == validate ]] && exit 0

index="${2:-${LSB_JOBINDEX:-}}"
[[ "$index" =~ ^[0-9]+$ ]] && (( index >= 1 && index <= 94 )) || usage
command -v genetribe >/dev/null || { echo 'GeneTribe executable not found in PATH' >&2; exit 1; }
query="${species[$((index-1))]}"
mkdir -p "$WORK_DIR/$query"
query_dir="$(cd "$WORK_DIR/$query" && pwd)"
cd "$query_dir"

link_input() {
  local sp="$1" ext="$2" target="$DATA_DIR/$sp.$ext" link="$query_dir/$sp.$ext"
  if [[ -L "$link" ]]; then
    [[ "$(readlink -f "$link")" == "$target" ]] || { echo "Conflicting link: $link" >&2; exit 1; }
  elif [[ -e "$link" ]]; then
    [[ "$link" -ef "$target" ]] || { echo "Conflicting file: $link" >&2; exit 1; }
  else
    ln -s "$target" "$link"
  fi
}

for ext in fa bed chrlist; do link_input "$query" "$ext"; done
for ref in "${REFS[@]}"; do
  [[ "${query,,}" == "${ref,,}" ]] && continue
  for ext in fa bed chrlist; do link_input "$ref" "$ext"; done
  outdir="ref_${ref}_output"
  if [[ -e "$outdir/.complete" ]]; then
    echo "Skipping completed comparison: $query versus $ref"
    continue
  fi
  [[ ! -e "$outdir" ]] || { echo "Incomplete output exists: $query_dir/$outdir; inspect before rerunning" >&2; exit 1; }
  [[ ! -e genetribe_output ]] || { echo "Unclaimed genetribe_output in $query_dir" >&2; exit 1; }
  for suffix in "${SUFFIXES[@]}"; do
    for base in "${query}_${ref}" "${ref}_${query}"; do
      [[ ! -e "$base.$suffix" ]] || { echo "Unclaimed GeneTribe result: $base.$suffix" >&2; exit 1; }
    done
  done
  mkdir "$outdir"
  echo "Running GeneTribe: reference=$ref query=$query"
  genetribe core -l "$ref" -f "$query" > "$outdir/runtime.log" 2>&1
  for suffix in "${SUFFIXES[@]}"; do
    for base in "${query}_${ref}" "${ref}_${query}"; do
      [[ ! -e "$base.$suffix" ]] || mv "$base.$suffix" "$outdir/"
    done
  done
  [[ ! -d genetribe_output ]] || mv genetribe_output "$outdir/"
  touch "$outdir/.complete"
done
echo "Completed query $index: $query"
