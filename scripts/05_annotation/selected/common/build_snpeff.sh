#!/usr/bin/env bash
#SBATCH --job-name=grass_snpeff_build
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
#SBATCH --output=logs/snpeff_build_%j.out
#SBATCH --error=logs/snpeff_build_%j.err
set -euo pipefail
: "${SNPEFF_IMAGE:?}" "${SNPEFF_CONFIG:?}" "${SNPEFF_DB:?}"
[[ -s "$SNPEFF_IMAGE" && -s "$SNPEFF_CONFIG" ]] || exit 1
# The historical jobs appended the same genome entry on every run. Add it once.
if ! grep -Eq "^[[:space:]]*${SNPEFF_DB}\.genome[[:space:]]*:" "$SNPEFF_CONFIG"; then
  printf '\n%s.genome : %s\n' "$SNPEFF_DB" "$SNPEFF_DB" >> "$SNPEFF_CONFIG"
fi
singularity exec "$SNPEFF_IMAGE" java -jar /opt/snpEff/snpEff.jar \
  build -c "$SNPEFF_CONFIG" -gff3 -v "$SNPEFF_DB"
singularity exec "$SNPEFF_IMAGE" java -jar /opt/snpEff/snpEff.jar \
  dump -c "$SNPEFF_CONFIG" "$SNPEFF_DB" > "${SNPEFF_DB}.dump.txt"
