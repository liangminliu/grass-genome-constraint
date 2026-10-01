#!/usr/bin/env bash
#SBATCH --job-name=atau_gerp_snp_intersect
#SBATCH --cpus-per-task=10
#SBATCH --array=1-4%4
#SBATCH --output=logs/atau_gerp_%A_%a.out
#SBATCH --error=logs/atau_gerp_%A_%a.err

# Four cumulative GERP thresholds: >0 is historical/exploratory; >2, >4,
# and >6 correspond to the manuscript's principal constraint classes.
set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit as an array job}"
VCF=${VCF:-141sp.filtered.bi-allelic.snp.vcf.gz}
GERP_BASE=${GERP_BASE:-/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Atau_msa/results/GERP_final_output}
IMAGE=${IMAGE:-$HOME/liuliangmin/4.may/software/Reseq_genek.sif}
thresholds=(merged_gt0 merged_gt2 merged_gt4 merged_gt6)
idx=$((SLURM_ARRAY_TASK_ID - 1))
[[ $idx -ge 0 && $idx -lt ${#thresholds[@]} ]] || exit 2
name=${thresholds[$idx]}
gerp="$GERP_BASE/${name}.bed"
[[ -s "$VCF" && -s "$gerp" && -s "$IMAGE" ]] || {
  echo "Missing VCF, GERP BED, or container image" >&2; exit 1;
}
for cmd in singularity bedtools bgzip; do command -v "$cmd" >/dev/null || exit 1; done
mkdir -p "$name" logs
clean=$(mktemp "$name/${name}.bed.XXXXXX")
trap 'rm -f "$clean"' EXIT
tail -n +2 "$gerp" > "$clean"  # original GERP BED includes a header

out_vcf="$name/${name}_snps.vcf.gz"
out_snp_bed="$name/${name}_snps.bed"
out_gerp_bed="$name/${name}_gerp_regions.bed"
singularity exec "$IMAGE" bcftools view "$VCF" |
  bedtools intersect -header -u -a stdin -b "$clean" |
  bgzip -c > "$out_vcf"
singularity exec "$IMAGE" bcftools index -f "$out_vcf"
singularity exec "$IMAGE" bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' "$out_vcf" > "$out_snp_bed"
bedtools intersect -u -a "$clean" -b "$out_snp_bed" > "$out_gerp_bed"
echo "Wrote $name SNP VCF/BED and intersecting GERP regions"
