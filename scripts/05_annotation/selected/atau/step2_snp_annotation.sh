#!/usr/bin/env bash
#SBATCH --job-name=atau_region_annotation
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/atau_region_%A_%a.out
#SBATCH --error=logs/atau_region_%A_%a.err

# Intersect the filtered SNP panel with one functional BED per array task.
# Submit with --array=1-N, where N is the number of lines in BED_LIST.
set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit as an array job}"
VCF=${VCF:-141sp.filtered.bi-allelic.snp.vcf.gz}
BED_LIST=${BED_LIST:-bed_list.txt}
IMAGE=${IMAGE:-$HOME/liuliangmin/4.may/software/Reseq_genek.sif}
[[ -s "$VCF" && -s "$BED_LIST" && -s "$IMAGE" ]] || {
  echo "Missing VCF, BED list, or container image" >&2; exit 1;
}
for cmd in singularity bedtools; do command -v "$cmd" >/dev/null || exit 1; done
bed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BED_LIST")
[[ -s "$bed" ]] || { echo "No readable BED at task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
name=$(basename "$bed" .bed)
mkdir -p output_bed output_vcf logs

singularity exec "$IMAGE" bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' "$VCF" |
  bedtools intersect -wa -wb -a stdin -b "$bed" > "output_bed/${name}.intersect.bed"
singularity exec "$IMAGE" bcftools view -R "$bed" -Oz \
  -o "output_vcf/${name}.vcf.gz" "$VCF"
singularity exec "$IMAGE" bcftools index -f "output_vcf/${name}.vcf.gz"
echo "Wrote region annotation for $name"
