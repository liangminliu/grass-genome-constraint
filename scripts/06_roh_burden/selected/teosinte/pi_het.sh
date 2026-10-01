#!/usr/bin/env bash
#SBATCH --job-name=zmay_pi_het
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/zmay_pi_het_%j.out
#SBATCH --error=logs/zmay_pi_het_%j.err
set -euo pipefail
VCF=${VCF:-all.filtered.bi-allelic.miss0.5.snp.vcf.gz}
OUTDIR=${OUTDIR:-vcftools_miss0.5_results}
PREFIX=${PREFIX:-all.filtered.bi-allelic.miss0.5}
export VCF OUTDIR PREFIX
exec bash "$(dirname "$0")/../common/calculate_pi_heterozygosity.sh"
