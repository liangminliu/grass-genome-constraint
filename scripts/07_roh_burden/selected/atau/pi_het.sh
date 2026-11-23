#!/usr/bin/env bash
#SBATCH --job-name=atau_pi_het
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/atau_pi_het_%j.out
#SBATCH --error=logs/atau_pi_het_%j.err
set -euo pipefail
VCF=${VCF:-141sp.filtered.bi-allelic.snp.vcf.gz}
OUTDIR=${OUTDIR:-vcftools_141sp_results}
PREFIX=${PREFIX:-141sp.filtered.bi-allelic}
export VCF OUTDIR PREFIX
exec bash "$(dirname "$0")/../common/calculate_pi_heterozygosity.sh"
