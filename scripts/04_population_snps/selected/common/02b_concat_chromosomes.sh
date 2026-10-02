#!/usr/bin/env bash
# Usage: 02b_concat_chromosomes.sh <chromosome_vcf_list.txt> <out.vcf.gz>
# The list contains one absolute path to a per-chromosome PASS SNP VCF per line,
# in reference chromosome order.
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Expected ordered VCF list and output VCF' >&2; exit 2; }
list="$1"; out="$2"
[[ -s "$list" ]] || { echo "Missing VCF list: $list" >&2; exit 1; }
[[ ! -e "$out" ]] || { echo "Existing output: $out" >&2; exit 1; }
while IFS= read -r path; do
  [[ -z "$path" ]] && continue
  [[ -s "$path" ]] || { echo "Missing chromosome VCF: $path" >&2; exit 1; }
done < "$list"
mkdir -p "$(dirname "$out")"
bcftools concat -f "$list" -Oz -o "$out"
bcftools index -f -t "$out"
