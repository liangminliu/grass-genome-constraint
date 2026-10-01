#!/usr/bin/env bash
#SBATCH --job-name=grass_sift_annotate
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=40
#SBATCH --mem=400G
#SBATCH --output=logs/sift_annotate_%j.out
#SBATCH --error=logs/sift_annotate_%j.err
set -euo pipefail
: "${SIFT4G_JAR:?}" "${DB_DIR:?}" "${VCF_IN:?}" "${OUT_DIR:?}"
[[ -s "$SIFT4G_JAR" && -d "$DB_DIR" && -s "$VCF_IN" ]] || exit 1
mkdir -p "$OUT_DIR/logs"
tmp=$(mktemp "$OUT_DIR/sift_input.XXXXXX.vcf")
trap 'rm -f "$tmp"' EXIT
if [[ "$VCF_IN" == *.gz ]]; then gzip -dc "$VCF_IN" > "$tmp"; else cat "$VCF_IN" > "$tmp"; fi
java -Xms"${JAVA_MEMORY:-120G}" -Xmx"${JAVA_MEMORY:-120G}" -jar "$SIFT4G_JAR" \
  -c -i "$tmp" -d "$DB_DIR" -r "$OUT_DIR" -t > "$OUT_DIR/logs/sift4g.log" 2>&1
grep -qi 'Completed successfully' "$OUT_DIR/logs/sift4g.log" || { echo 'Check SIFT4G log' >&2; exit 1; }
