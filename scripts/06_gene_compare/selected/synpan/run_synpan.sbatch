#!/usr/bin/env bash
#SBATCH --job-name=SynPan_grass
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --output=logs/synpan_%j.out
#SBATCH --error=logs/synpan_%j.err
set -euo pipefail
: "${INPUT_DIR:?Set INPUT_DIR to the directory containing species .prot and .bed files}"
: "${WORK_DIR:?Set WORK_DIR to a dedicated writable SynPan run directory}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[[ "$INPUT_DIR" = /* && "$WORK_DIR" = /* ]] || { echo 'Use absolute input and work paths' >&2; exit 1; }
[[ "$INPUT_DIR" != "$WORK_DIR" ]] || { echo 'Keep input and work directories separate' >&2; exit 1; }
mkdir -p "$WORK_DIR"
[[ ! -e "$WORK_DIR/grass.list.SG" && ! -e "$WORK_DIR/Atau.prot_raw" ]] || {
  echo 'Use a fresh WORK_DIR for each SynPan run' >&2; exit 1;
}
cp "$script_dir/grass.list" "$WORK_DIR/grass.list"
while IFS= read -r species; do
  [[ "$species" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 1
  for ext in prot bed; do
    source_file=$INPUT_DIR/$species.$ext
    [[ -s "$source_file" ]] || { echo "Missing input: $source_file" >&2; exit 1; }
    target=$WORK_DIR/$species.$ext
    if [[ ! -e "$target" && ! -L "$target" ]]; then ln -s "$source_file" "$target"; fi
  done
done < "$script_dir/grass.list"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-40}
cd "$WORK_DIR"
perl "$script_dir/../synpan_build.pl" grass.list
[[ -s grass.list.SG && -s grass.list.SG.pan && -s grass.list.sg.stat ]] || {
  echo 'SynPan outputs incomplete' >&2; exit 1;
}
echo "SynPan finished at $(date)"
