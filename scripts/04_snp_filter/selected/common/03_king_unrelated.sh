#!/usr/bin/env bash
# Infer kinship before selecting the final bamboo or teosinte sample set.
# Usage: KING_DEGREE=<recorded_degree> 03_king_unrelated.sh <hard-filtered.vcf.gz> <out_prefix>
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Expected hard-filtered VCF and output prefix' >&2; exit 2; }
[[ "${KING_DEGREE:-}" =~ ^[1-9][0-9]*$ ]] || {
  echo 'Set KING_DEGREE to the value documented in the production run' >&2; exit 2;
}
input="$1"; out="$2"
[[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
for tool in plink king; do command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }; done
[[ ! -e "$out.bed" ]] || { echo "Existing PLINK output: $out.bed" >&2; exit 1; }
mkdir -p "$(dirname "$out")"
plink --vcf "$input" --make-bed --allow-extra-chr --out "$out"
king -b "$out.bed" --kinship --prefix "$out.kinship" --cpus "${THREADS:-8}"
king -b "$out.bed" --unrelated --degree "$KING_DEGREE" \
  --prefix "$out.unrelated" --cpus "${THREADS:-8}"
echo 'Inspect KING output and create the approved one-column final keep list before final filtering.'
