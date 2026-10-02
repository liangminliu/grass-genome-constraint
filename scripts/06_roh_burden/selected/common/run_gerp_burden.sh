#!/usr/bin/env bash
# Usage: run_gerp_burden.sh <ped|teosinte|atau> <froh.tsv> <roh.tsv> <output>
# Set the species-specific VCF variables described below before running.
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Expected species FROH ROH output' >&2; exit 2; }
species=$1; froh=$2; roh=$3; output=$4
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
for f in "$froh" "$roh"; do [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }; done

case "$species" in
  ped)
    : "${PEDC_GERP4:?Set PEDC_GERP4}"
    : "${PEDD_GERP4:?Set PEDD_GERP4}"
    : "${GERP_LT0:?Set GERP_LT0}"
    : "${ANCESTRAL_BED:?Set ANCESTRAL_BED}"
    : "${CHROM_MAP:?Set CHROM_MAP}"
    for f in "$PEDC_GERP4" "$PEDD_GERP4" "$GERP_LT0" "$ANCESTRAL_BED" "$CHROM_MAP"; do [[ -s "$f" ]] || exit 1; done
    mkdir -p "$output"
    python3 "$root/ped/compute_gerp_burden.py" \
      --vcf "GERP4:PedC:$PEDC_GERP4" --vcf "GERP4:PedD:$PEDD_GERP4" \
      --vcf "GERP_lt0:Ped_lt0:$GERP_LT0" --ancestral-bed "$ANCESTRAL_BED" \
      --froh "$froh" --roh "$roh" --chrom-map "$CHROM_MAP" \
      --chrom-normalize auto --outdir "$output"
    ;;
  teosinte)
    : "${GERP4_VCF:?Set GERP4_VCF}"
    : "${GERP_LT0:?Set GERP_LT0}"
    for f in "$GERP4_VCF" "$GERP_LT0"; do [[ -s "$f" ]] || exit 1; done
    mkdir -p "$output"
    python3 "$root/teosinte/compute_gerp_burden.py" \
      --vcf "GERP4:Teosinte_gt4:$GERP4_VCF" \
      --vcf "GERP_lt0:Teosinte_lt0:$GERP_LT0" \
      --froh "$froh" --roh "$roh" --chrom-normalize auto --outdir "$output"
    ;;
  atau)
    : "${GERP4_VCF:?Set GERP4_VCF}"
    : "${GERP_LT0:?Set GERP_LT0}"
    : "${LAUTO_FILE:?Set LAUTO_FILE}"
    for f in "$GERP4_VCF" "$GERP_LT0" "$LAUTO_FILE"; do [[ -s "$f" ]] || exit 1; done
    mkdir -p "$(dirname "$output")"
    python3 "$root/atau/compute_gerp_burden.py" \
      --gerp4-vcf "$GERP4_VCF" --gerpneg-vcf "$GERP_LT0" \
      --froh "$froh" --roh "$roh" --lauto "$LAUTO_FILE" \
      --out-prefix "$output" --exclude-samples "${EXCLUDE_SAMPLES:-BW_01192}" --ploidy 2
    ;;
  *) echo "Unknown species: $species" >&2; exit 2 ;;
esac
