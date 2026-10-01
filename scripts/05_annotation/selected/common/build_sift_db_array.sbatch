#!/usr/bin/env bash
#SBATCH --job-name=grass_sift_db
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=40
#SBATCH --mem=300G
#SBATCH --output=logs/sift_db_%A_%a.out
#SBATCH --error=logs/sift_db_%A_%a.err
set -euo pipefail
for name in GENOME_DIR ORG ORG_VERSION SIFT_DIR SIFT4G_BIN PROTEIN_DB SLURM_ARRAY_TASK_ID; do
  [[ -n "${!name:-}" ]] || { echo "Missing $name" >&2; exit 1; }
done
list=${CHR_LIST:-$GENOME_DIR/configs/chr_allow.list}
fasta=${FASTA_IN:-$GENOME_DIR/chr-src/$ORG.fa}
gtf=${GTF_IN:-$GENOME_DIR/gene-annotation-src/$ORG.gtf.gz}
for file in "$list" "$fasta" "$fasta.fai" "$gtf" "$SIFT_DIR/make-SIFT-db-all.pl" "$SIFT4G_BIN" "$PROTEIN_DB"; do
  [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
mapfile -t chrs < <(awk 'NF && $1 !~ /^#/ {print $1}' "$list")
idx=$((SLURM_ARRAY_TASK_ID - 1))
(( idx >= 0 && idx < ${#chrs[@]} )) || exit 1
chr=${chrs[$idx]}
[[ "$chr" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 1
work=$GENOME_DIR/build_$chr
mkdir -p "$work"/{chr-src,gene-annotation-src,logs,configs}
samtools faidx "$fasta" "$chr" | gzip -c > "$work/chr-src/$chr.fa.gz"
zcat "$gtf" | awk -v c="$chr" '$1 == c' | gzip -c > "$work/gene-annotation-src/$ORG.$chr.gtf.gz"
records=$(gzip -dc "$work/gene-annotation-src/$ORG.$chr.gtf.gz" | awk 'END {print NR+0}')
(( records > 0 )) || { echo "No annotation records for $chr" >&2; exit 1; }
config=$work/configs/config.$chr.txt
cat > "$config" <<EOF
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=1
MITO_GENETIC_CODE_TABLENAME=Plant_mito_uses_Standard
PLASTID_GENETIC_CODE_TABLE=11
PLASTID_GENETIC_CODE_TABLENAME=Plant_Plastid
PARENT_DIR=$work
ORG=$ORG
ORG_VERSION=$ORG_VERSION
SIFT4G_PATH=$SIFT4G_BIN
PROTEIN_DB=$PROTEIN_DB
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord
EOF
cd "$SIFT_DIR"
perl make-SIFT-db-all.pl -config "$config" > "$work/logs/run_$chr.log" 2>&1
grep -q 'All done' "$work/logs/run_$chr.log" || exit 1
touch "$work/.SUCCESS"
