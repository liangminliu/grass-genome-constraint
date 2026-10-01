#!/usr/bin/env bash
# Paths follow the recovered Atau jobs and author-supplied SnpEff commands.
SIFT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Atau_msa/pop/04_SIFT/scripts_to_build_SIFT_db
GENOME_DIR=$SIFT_DIR/test_files/Atau
ORG=Atau
ORG_VERSION=Atau.V1
CHR_LIST=$GENOME_DIR/configs/chr_allow.list
FASTA_IN=$GENOME_DIR/chr-src/Atau.fa
GTF_IN=$GENOME_DIR/gene-annotation-src/Atau.gtf.gz
SIFT4G_BIN=/gpfs/hpc/home/kib_unit/lidezhu/miniforge3/envs/sift4g/bin/sift4g
PROTEIN_DB=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/4.may/step3_genetic_load/2.database/scripts_to_build_SIFT_db/databases/uniref90.fasta
DB_DIR=$GENOME_DIR/Atau_full
SIFT4G_JAR=/gpfs/hpc/home/kib_unit/lidezhu/zf/sift4g/SIFT4G_Annotator.jar
SIFT_VCF=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Atau_msa/pop/04_SIFT/SIFI_annotator/141sp.filtered.bi-allelic.snp.vcf.gz
SIFT_OUT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/2.msa/Atau_msa/pop/04_SIFT/SIFI_annotator/SIFT4G_results
SNPEFF_IMAGE=$HOME/liuliangmin/4.may/software/Reseq_genek.sif
SNPEFF_CONFIG=./snpEff.config
SNPEFF_DB=Atau
SNPEFF_ALL_VCF=141sp.filtered.bi-allelic.snp.vcf.gz
SNPEFF_ALL_OUT=141sp.filtered.bi-allelic.snp.ann.vcf
GERP_DIR=GERP_del
GERP_LABELS=merged
SIFT_DELETERIOUS_VCF=SIFT_del/Atau.sift_lt_all.vcf.gz
SIFT_DELETERIOUS_OUT=SIFT_del/Atau.sift_lt_all.ann.vcf
