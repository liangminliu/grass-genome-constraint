#!/usr/bin/env bash
#SBATCH --job-name=atau_vcf_filter
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/atau_vcf_filter.%j.out
#SBATCH --error=logs/atau_vcf_filter.%j.err

# Author-supplied A. tauschii 141-sample pre-final VCF filter.
# The input must already contain SNPs; the original command does not remove indels.
set -euo pipefail
input=${INPUT_VCF:-extracted_141sp_Atau.vcf.gz}
output=${OUTPUT_VCF:-141sp.filtered.bi-allelic.snp.vcf.gz}
[[ -s "$input" ]] || { echo "Missing input VCF: $input" >&2; exit 1; }
for cmd in vcftools bgzip tabix; do
  command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 1; }
done
mkdir -p logs
vcftools --gzvcf "$input" \
  --minDP 1 --maxDP 100 --minGQ 10 --minQ 30 \
  --max-missing 0.5 --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all --stdout | bgzip -c > "$output"
tabix -f -p vcf "$output"
echo "Wrote $output"
