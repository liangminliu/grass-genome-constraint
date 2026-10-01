#!/usr/bin/env bash
#SBATCH --job-name=ped_pi_het
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/ped_pi_het_%j.out
#SBATCH --error=logs/ped_pi_het_%j.err
set -euo pipefail
VCF=${VCF:-63sp.filtered.mainChr.vcf.gz}
OUTDIR=${OUTDIR:-vcftools_63sp_results}
PREFIX=${PREFIX:-63sp.filtered.mainChr}
export VCF OUTDIR PREFIX
exec bash "$(dirname "$0")/../common/calculate_pi_heterozygosity.sh"
