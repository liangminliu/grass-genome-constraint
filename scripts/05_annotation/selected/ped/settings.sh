#!/usr/bin/env bash
# Paths follow the recovered Ped/PED jobs. Supply the referenced inputs on the cluster.
SIFT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/5.deleterious/06_SIFT/scripts_to_build_SIFT_db
GENOME_DIR=$SIFT_DIR/test_files/PED
ORG=PED
ORG_VERSION=PED.V2
CHR_LIST=$GENOME_DIR/configs/chr_allow.list
FASTA_IN=$GENOME_DIR/chr-src/ped.fa
GTF_IN=$GENOME_DIR/gene-annotation-src/PED.chrOnly.gtf.gz
SIFT4G_BIN=/gpfs/hpc/home/kib_unit/lidezhu/miniforge3/envs/sift4g/bin/sift4g
PROTEIN_DB=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/4.may/step3_genetic_load/2.database/scripts_to_build_SIFT_db/databases/uniref90.fasta
DB_DIR=$GENOME_DIR/PED.V2.full
SIFT4G_JAR=/gpfs/hpc/home/kib_unit/lidezhu/zf/sift4g/SIFT4G_Annotator.jar
SIFT_VCF=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/5.deleterious/06_SIFT/SIFI_annotator/63sp.filtered.mainChr.vcf.gz
SIFT_OUT_DIR=/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/5.deleterious/06_SIFT/SIFI_annotator/SIFT4G_results/ALL_annot
SNPEFF_IMAGE=$HOME/liuliangmin/4.may/software/Reseq_genek.sif
SNPEFF_CONFIG=./snpEff.config
SNPEFF_DB=ped
SNPEFF_ALL_VCF=71sp.filtered.bi-allelic.snp.rename.vcf.gz
SNPEFF_ALL_OUT=71sp.filtered.snp.ann.vcf
GERP_DIR=GERP_del
GERP_LABELS='pedC pedD'
