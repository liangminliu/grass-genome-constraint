#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Usage: workflow.sh settings.sh {sift-build|sift-merge|sift-annotate|snpeff-build|snpeff-all|snpeff-gerp|snpeff-sift}' >&2; exit 2; }
settings=$(realpath "$1")
action=$2
# Species settings are local path declarations maintained with this release.
source "$settings"
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p logs
case "$action" in
  sift-build)
    mkdir -p logs
    count=$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "${CHR_LIST:-$GENOME_DIR/configs/chr_allow.list}")
    (( count > 0 )) || exit 1
    export GENOME_DIR ORG ORG_VERSION SIFT_DIR SIFT4G_BIN PROTEIN_DB CHR_LIST FASTA_IN GTF_IN
    sbatch --export=ALL --array="1-$count" "$here/build_sift_db_array.sbatch"
    ;;
  sift-merge) export GENOME_DIR ORG_VERSION DB_DIR CHR_LIST; bash "$here/merge_sift_db.sh" ;;
  sift-annotate)
    export SIFT4G_JAR DB_DIR
    VCF_IN=${SIFT_VCF:?}; OUT_DIR=${SIFT_OUT_DIR:?}; export VCF_IN OUT_DIR
    sbatch --export=ALL "$here/annotate_sift4g.sh"
    ;;
  snpeff-build) export SNPEFF_IMAGE SNPEFF_CONFIG SNPEFF_DB; sbatch --export=ALL "$here/build_snpeff.sh" ;;
  snpeff-all)
    VCF_IN=${SNPEFF_ALL_VCF:?}; VCF_OUT=${SNPEFF_ALL_OUT:?}; export VCF_IN VCF_OUT SNPEFF_IMAGE SNPEFF_CONFIG SNPEFF_DB
    sbatch --export=ALL "$here/annotate_snpeff.sh"
    ;;
  snpeff-gerp)
    export SNPEFF_IMAGE SNPEFF_CONFIG SNPEFF_DB
    for label in ${GERP_LABELS:?}; do
      VCF_IN=${GERP_DIR:?}/${label}_gt4_snps.vcf.gz
      VCF_OUT=${GERP_DIR}/${label}_gt4.snp.ann.vcf
      export VCF_IN VCF_OUT
      sbatch --export=ALL "$here/annotate_snpeff.sh"
    done
    ;;
  snpeff-sift)
    VCF_IN=${SIFT_DELETERIOUS_VCF:?}; VCF_OUT=${SIFT_DELETERIOUS_OUT:?}; export VCF_IN VCF_OUT SNPEFF_IMAGE SNPEFF_CONFIG SNPEFF_DB
    sbatch --export=ALL "$here/annotate_snpeff.sh"
    ;;
  *) echo "Unknown action: $action" >&2; exit 2 ;;
esac
