#!/usr/bin/env bash
# Methods-based reconstruction of the executed bamboo/teosinte read-to-gVCF stage.
# Usage: 01_reads_to_gvcf.sh <sample> <R1.fastq.gz> <R2.fastq.gz> <reference.fa> <outdir>
set -euo pipefail
[[ $# -eq 5 ]] || { echo 'Expected sample R1 R2 reference outdir' >&2; exit 2; }
sample="$1"; r1="$2"; r2="$3"; ref="$4"; out="$5"
for f in "$r1" "$r2" "$ref"; do [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }; done
for tool in fastp bwa samtools picard gatk; do command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }; done
mkdir -p "$out"
[[ ! -e "$out/$sample.g.vcf.gz" ]] || { echo "Existing gVCF: $sample" >&2; exit 1; }
fastp -i "$r1" -I "$r2" -o "$out/$sample.R1.clean.fq.gz" -O "$out/$sample.R2.clean.fq.gz" \
  --detect_adapter_for_pe --trim_poly_g --qualified_quality_phred 20 --unqualified_percent_limit 40 \
  --length_required 70 --json "$out/$sample.fastp.json" --html "$out/$sample.fastp.html"
bwa mem -t "${THREADS:-8}" -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$ref" \
  "$out/$sample.R1.clean.fq.gz" "$out/$sample.R2.clean.fq.gz" \
  | samtools sort -@ "${THREADS:-8}" -o "$out/$sample.sorted.bam" -
samtools view -@ "${THREADS:-8}" -b -q 20 -o "$out/$sample.mapq20.bam" "$out/$sample.sorted.bam"
picard MarkDuplicates I="$out/$sample.mapq20.bam" O="$out/$sample.dedup.bam" \
  M="$out/$sample.duplicate_metrics.txt" REMOVE_DUPLICATES=true CREATE_INDEX=true
gatk HaplotypeCaller -R "$ref" -I "$out/$sample.dedup.bam" -O "$out/$sample.g.vcf.gz" -ERC GVCF
