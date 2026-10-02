#!/usr/bin/env bash
# Usage: 02_joint_call_chromosome.sh <chrom> <reference.fa> <sample_map.tsv> <outdir>
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Expected chromosome reference sample_map outdir' >&2; exit 2; }
chrom="$1"; ref="$2"; sample_map="$3"; out="$4"
[[ "$chrom" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Unsafe chromosome name' >&2; exit 1; }
for f in "$ref" "$sample_map"; do [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }; done
command -v gatk >/dev/null || { echo 'Missing GATK' >&2; exit 1; }
mkdir -p "$out"
[[ ! -e "$out/$chrom.pass.snps.vcf.gz" && ! -e "$out/$chrom.gdb" ]] || { echo "Existing output for $chrom" >&2; exit 1; }
gatk GenomicsDBImport --sample-name-map "$sample_map" --genomicsdb-workspace-path "$out/$chrom.gdb" \
  --intervals "$chrom" --reader-threads "${THREADS:-8}"
gatk GenotypeGVCFs -R "$ref" -V "gendb://$out/$chrom.gdb" -O "$out/$chrom.raw.vcf.gz"
gatk SelectVariants -R "$ref" -V "$out/$chrom.raw.vcf.gz" --select-type-to-include SNP \
  -O "$out/$chrom.snps.vcf.gz"
gatk VariantFiltration -R "$ref" -V "$out/$chrom.snps.vcf.gz" \
  --filter-name MethodsHardFilter \
  --filter-expression 'QD < 2.0 || QUAL < 30.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
  -O "$out/$chrom.filtered.snps.vcf.gz"
bcftools view -f PASS -Oz -o "$out/$chrom.pass.snps.vcf.gz" "$out/$chrom.filtered.snps.vcf.gz"
bcftools index -f -t "$out/$chrom.pass.snps.vcf.gz"
