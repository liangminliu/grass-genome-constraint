#!/usr/bin/env bash
# Release entry point for an existing msa_pipeline checkout, not a historical job.
# Usage: run_msa_with_existing_tree.sh WORKDIR CONFIG GUIDE_TREE [all|roast|conservation]
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 WORKDIR CONFIG GUIDE_TREE [all|roast|conservation]" >&2
  exit 2
fi
workdir=$(cd "$1" && pwd)
config=$(realpath "$2")
tree=$(realpath "$3")
stage=${4:-roast}
case "$stage" in all|roast|conservation) ;; *) echo "Invalid stage: $stage" >&2; exit 2;; esac
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[[ -s "$config" && -s "$tree" && -f "$workdir/workflow/Snakefile" ]] || {
  echo "Missing config, tree, or workflow/Snakefile" >&2; exit 1;
}
command -v snakemake >/dev/null || { echo "snakemake not found" >&2; exit 1; }
runtime_config="$workdir/config/release_existing_tree.yaml"
python3 "$script_dir/prepare_existing_tree_config.py" \
  --config "$config" --tree "$tree" --output "$runtime_config"
mkdir -p "$workdir/logs"
cd "$workdir"

run_step() {
  local target=$1 jobs=$2 cluster=$3
  echo "Starting $target with guide tree $tree"
  snakemake --snakefile workflow/Snakefile --configfile "$runtime_config" \
    --rerun-incomplete --jobs "$jobs" --cluster "$cluster" "$target" \
    > "logs/release_${target}.log" 2>&1 || {
      echo "Failed $target; see logs/release_${target}.log" >&2; exit 1;
    }
}

if [[ "$stage" == all ]]; then
  run_step align 50 'sbatch -p parallel -N 1 --ntasks=1 --cpus-per-task=4'
fi
if [[ "$stage" == all || "$stage" == roast ]]; then
  # Rebuild topology from supplied Newick even if an older topology.tre exists.
  snakemake --snakefile workflow/Snakefile --configfile "$runtime_config" \
    --cores 1 --forcerun write_tree results/tree/topology.tre \
    > logs/release_tree.log 2>&1 || {
      echo "Failed guide-tree installation; see logs/release_tree.log" >&2; exit 1;
    }
  run_step roast 1 'sbatch -p FAT1 -N 1 --ntasks=1 --cpus-per-task=60'
fi
if [[ "$stage" == all || "$stage" == conservation ]]; then
  run_step call_conservation 10 'sbatch -p parallel -N 1 --ntasks=1 --cpus-per-task=20'
fi
