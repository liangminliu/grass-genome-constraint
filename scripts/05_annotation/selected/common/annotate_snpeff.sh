#!/usr/bin/env bash
#SBATCH --job-name=grass_snpeff_annotate
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
#SBATCH --output=logs/snpeff_annotate_%j.out
#SBATCH --error=logs/snpeff_annotate_%j.err
set -euo pipefail
: "${SNPEFF_IMAGE:?}" "${SNPEFF_CONFIG:?}" "${SNPEFF_DB:?}" "${VCF_IN:?}" "${VCF_OUT:?}"
[[ -s "$SNPEFF_IMAGE" && -s "$SNPEFF_CONFIG" && -s "$VCF_IN" ]] || exit 1
mkdir -p "$(dirname "$VCF_OUT")"
prefix=${STATS_PREFIX:-${VCF_OUT%.vcf}}
singularity exec "$SNPEFF_IMAGE" java -Xmx"${JAVA_MEMORY:-200G}" -jar /opt/snpEff/snpEff.jar \
  -c "$SNPEFF_CONFIG" -ud 5000 \
  -csvStats "${prefix}.csv" -htmlStats "${prefix}.html" \
  -o vcf "$SNPEFF_DB" "$VCF_IN" > "$VCF_OUT"
bgzip -f "$VCF_OUT"
tabix -f -p vcf "$VCF_OUT.gz"
echo "Annotated VCF: $VCF_OUT.gz"
