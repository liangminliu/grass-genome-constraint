#!/usr/bin/env bash
# Usage: 03_final_population_filter.sh <ped|zmay> <hard-filtered.vcf.gz> <keep_samples.txt> <out_prefix>
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Expected species input_vcf keep_samples out_prefix' >&2; exit 2; }
species="$1"; input="$2"; keep="$3"; out="$4"
case "$species" in ped) expected=63 ;; zmay) expected=90 ;; *) echo 'Species must be ped or zmay' >&2; exit 2 ;; esac
for f in "$input" "$keep"; do [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }; done
for tool in bcftools vcftools bgzip tabix; do command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }; done
[[ ! -e "$out.vcf.gz" ]] || { echo "Existing output: $out.vcf.gz" >&2; exit 1; }
[[ $(awk 'NF && !/^#/ {print $1}' "$keep" | sort -u | wc -l) -eq "$expected" ]] || {
  echo "Expected $expected distinct sample IDs in $keep" >&2; exit 1;
}
mkdir -p "$(dirname "$out")"
bcftools view -S "$keep" -m2 -M2 -v snps -Oz -o "$out.subset.vcf.gz" "$input"
bcftools index -f -t "$out.subset.vcf.gz"
[[ $(bcftools query -l "$out.subset.vcf.gz" | wc -l) -eq "$expected" ]] || {
  echo 'VCF sample count does not match the final sample list' >&2; exit 1;
}
# VCFtools genotype filters correspond to Methods DP 1–100 and GQ >=10.
# The 50% site call-rate criterion is applied after genotype masking.
vcftools --gzvcf "$out.subset.vcf.gz" --minDP 1 --maxDP 100 --minGQ 10 \
  --minQ 30 --max-missing 0.5 --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all --stdout | bgzip -c > "$out.vcf.gz"
tabix -f -p vcf "$out.vcf.gz"
