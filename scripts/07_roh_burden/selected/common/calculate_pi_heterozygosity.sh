#!/usr/bin/env bash
set -euo pipefail
: "${VCF:?Set VCF to the filtered SNP VCF}"
: "${OUTDIR:?Set OUTDIR for vcftools results}"
: "${PREFIX:?Set PREFIX for output filenames}"
[[ -s "$VCF" ]] || { echo "Missing VCF: $VCF" >&2; exit 1; }
command -v vcftools >/dev/null || exit 1
mkdir -p "$OUTDIR"
base=$OUTDIR/$PREFIX
if [[ "$VCF" == *.gz ]]; then input=(--gzvcf "$VCF"); else input=(--vcf "$VCF"); fi
vcftools "${input[@]}" --site-pi --out "$base"
vcftools "${input[@]}" --het --out "$base"
awk 'BEGIN{OFS="\t"; print "INDV","Heterozygosity"}
     NR>1 {if ($4>0) printf "%s\t%.6f\n",$1,($4-$2)/$4; else print $1,"NA"}' \
  "$base.het" > "$base.heterozygosity.tsv"
awk 'NR>1 {sum+=$3; n++} END {if (n) printf "%.8f\n",sum/n; else print "NA"}' \
  "$base.sites.pi" > "$base.mean_pi.txt"
echo "Saved $base.sites.pi, $base.het, $base.heterozygosity.tsv and $base.mean_pi.txt"
