# Population SNP filtering

`selected/atau/step1_vcf_filter.sh` implements the supplied A. tauschii filtered biallelic-SNP branch. `selected/atau/extract_retained_samples.sh` applies `metadata/atau_retained_samples.txt`, checks every name against the VCF header and writes an indexed retained-sample VCF. The list has exactly 141 distinct names, including `BW_01192`. The separate ROH workflow removes `BW_01192` for 140 individuals.

Final bamboo and teosinte mapping, joint genotyping, KING and sample-filter launchers are pending production records, so candidate jobs are withheld.
