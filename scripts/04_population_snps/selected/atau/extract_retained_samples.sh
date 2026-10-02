#!/usr/bin/env bash
set -euo pipefail
: "${INPUT_VCF:?Set INPUT_VCF to the filtered A. tauschii VCF}"
: "${OUTPUT_VCF:?Set OUTPUT_VCF to the retained-sample VCF.gz}"
: "${SAMPLES:?Set SAMPLES to the retained-sample ID list}"
samples="$SAMPLES"
[[ -s "$INPUT_VCF" && -s "$samples" ]] || exit 1
[[ "$OUTPUT_VCF" == *.vcf.gz ]] || { echo 'OUTPUT_VCF must end in .vcf.gz' >&2; exit 1; }
[[ $(sort "$samples" | uniq -d | wc -l) -eq 0 ]] || { echo 'Duplicate sample names' >&2; exit 1; }
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
bcftools query -l "$INPUT_VCF" > "$tmp"
missing=$(comm -23 <(sort "$samples") <(sort "$tmp"))
[[ -z "$missing" ]] || { printf 'Samples absent from VCF:\n%s\n' "$missing" >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT_VCF")"
bcftools view -S "$samples" -Oz -o "$OUTPUT_VCF" "$INPUT_VCF"
bcftools index -t -f "$OUTPUT_VCF"
echo "Retained $(wc -l < "$samples") samples in $OUTPUT_VCF"

