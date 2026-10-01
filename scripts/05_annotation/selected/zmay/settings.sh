#!/usr/bin/env bash
# Paths follow the recovered Zmay jobs.
SIFT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/5.deleterious/06_SIFT/scripts_to_build_SIFT_db
GENOME_DIR=$SIFT_DIR/test_files/Zmay
ORG=Zmay
ORG_VERSION=Zmay.V5
CHR_LIST=$GENOME_DIR/configs/chr_allow.list
FASTA_IN=$GENOME_DIR/chr-src/Zmay.fa
GTF_IN=$GENOME_DIR/gene-annotation-src/Zmay.gtf.gz
SIFT4G_BIN=/gpfs/hpc/home/kib_unit/lidezhu/miniforge3/envs/sift4g/bin/sift4g
PROTEIN_DB=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/4.may/step3_genetic_load/2.database/scripts_to_build_SIFT_db/databases/uniref90.fasta
DB_DIR=$GENOME_DIR/Zmay.V5.full
SIFT4G_JAR=/gpfs/hpc/home/kib_unit/lidezhu/zf/sift4g/SIFT4G_Annotator.jar
SIFT_VCF=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Zmay/Zmay_v5_results/GERP_output/03.deletrious/04.SIFT/all.filtered.bi-allelic.miss0.5.snp.vcf.gz
SIFT_OUT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Zmay/Zmay_v5_results/GERP_output/03.deletrious/04.SIFT/SIFT4G_results
SNPEFF_IMAGE=$HOME/liuliangmin/4.may/software/Reseq_genek.sif
SNPEFF_CONFIG=./snpEff.config
SNPEFF_DB=Zmay
SNPEFF_ALL_VCF=all.filtered.bi-allelic.miss0.5.snp.vcf.gz
SNPEFF_ALL_OUT=all.filtered.bi-allelic.miss0.5.snp.ann.vcf
GERP_DIR=GERP_del
GERP_LABELS=merged
SIFT_DELETERIOUS_VCF=SIFT_del/Zmay.sift_lt_all.vcf.gz
SIFT_DELETERIOUS_OUT=SIFT_del/Zmay.sift_lt_all.ann.vcf
