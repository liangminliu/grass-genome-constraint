#!/usr/bin/env bash
#SBATCH --job-name=maf_coverage
#SBATCH --output=logs/coverage_%A_%a.out
#SBATCH --error=logs/coverage_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# Reusable for Ped, Zmay, and Atau reference-coordinate alignments.
# Submit with --array=1-N, where N is the number of query IDs in SPECIES_LIST.
set -euo pipefail
: "${REF_SPECIES:?Set the reference sequence prefix (for example ped, Zmay, or Atau)}"
: "${SLURM_ARRAY_TASK_ID:?Submit this job as a SLURM array}"
MAF_FILE=${MAF_FILE:-roast.maf}
SPECIES_LIST=${SPECIES_LIST:-94sp.pop.list}
OUT_DIR=${OUT_DIR:-mafPairCoverage_results}
[[ -s "$MAF_FILE" && -s "$SPECIES_LIST" ]] || { echo "Missing MAF or species list" >&2; exit 1; }
command -v mafPairCoverage >/dev/null || { echo "mafPairCoverage is required" >&2; exit 1; }

target=$(awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print $1}' "$SPECIES_LIST")
[[ -n "$target" ]] || { echo "No query at row $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
[[ "$target" != "$REF_SPECIES" ]] || { echo "Reference appears in query list: $target" >&2; exit 1; }
[[ "$target" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Invalid query ID: $target" >&2; exit 1; }
mkdir -p "$OUT_DIR" logs
echo "Reference=$REF_SPECIES query=$target MAF=$MAF_FILE"
mafPairCoverage \
  --maf "$MAF_FILE" \
  --seq1 "${REF_SPECIES}*" \
  --seq2 "${target}*" \
  --bin_length 5000 \
  > "$OUT_DIR/${target}.tsv" \
  2> "logs/${target}.coverage.log"
[[ -s "$OUT_DIR/${target}.tsv" ]] || { echo "Empty coverage result for $target" >&2; exit 1; }
