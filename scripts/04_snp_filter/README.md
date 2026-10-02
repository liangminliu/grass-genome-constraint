# Population SNP filtering

`selected/atau/step1_vcf_filter.sh` implements the supplied A. tauschii filtered biallelic-SNP branch. `selected/atau/extract_retained_samples.sh` applies `metadata/atau_retained_samples.txt`, checks every name against the VCF header and writes an indexed retained-sample VCF. The list has exactly 141 distinct names, including `BW_01192`. The separate ROH workflow removes `BW_01192` for 140 individuals.

The exact bamboo/teosinte production job logs, final unrelated-sample lists and VCF hashes remain pending.

# Bamboo and teosinte SNP filtering

`selected/common/` contains a Methods-based, parameterized reconstruction of the bamboo and teosinte pipeline:

1. `01_reads_to_gvcf.sh`: fastp, BWA-MEM, MAPQ ≥20, duplicate removal and per-sample GATK HaplotypeCaller gVCF.
2. `02_joint_call_chromosome.sh` and `02b_concat_chromosomes.sh`: per-chromosome GenomicsDBImport, joint genotyping, SNP selection and the documented GATK hard filter, followed by concatenation in reference order.
3. `03_king_unrelated.sh`: convert the hard-filtered VCF for KING `--kinship`/`--unrelated`; set `KING_DEGREE` to the recorded production value and inspect the resulting unrelated set.
4. `04_final_population_filter.sh`: apply the approved unrelated-sample keep list (63 bamboo or 90 teosinte), retain biallelic SNPs, mask genotypes outside DP 1–100 or GQ <10, and retain sites with QUAL ≥30 and call rate ≥50%.

The scripts accept input paths explicitly and refuse to overwrite final outputs. The reference FASTA must be indexed for BWA, SAMtools and GATK; provide a two-column GATK sample-name map and a chromosome list matching the reference. The `ped` and `zmay` keep files must contain the exact final 63 and 90 IDs, respectively. KING relatedness inference precedes the final keep list: the manuscript specifies KING `--kinship` and `--unrelated`, but the exact production output and choice of unrelated set still need to be supplied. Consequently, these newly added upstream scripts are a **Methods reconstruction**, whereas the existing A. tauschii scripts are source-derived. Verify tool versions, VCF hashes and final sample identities against the production records before treating a rerun as the original analysis.
