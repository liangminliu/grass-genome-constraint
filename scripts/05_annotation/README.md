# Deleterious variant annotation

`selected/common/workflow.sh` runs the released three-species SIFT 4G and SnpEff build/annotation stages. Use one of `selected/{ped,zmay,atau}/settings.sh` from the species working directory and read `selected/common/README.md` for the commands. The workflow includes all-SNP and GERP>4 SnpEff outputs; bamboo C/D-specific SIFT classes are outside the shared wrapper and need the final class VCFs.

The A. tauschii `selected/atau/step2_snp_annotation.sh`, `step3_intersect_snps_gerp.sh` and `step4_functional_intersect.sh` supply the author-provided SNP-to-region and SNP-to-GERP route. The recorded SIFT/SnpEff annotation VCF name has 141 samples. Output/database hashes and the final SIFT <0.05 extraction command are requested in the root missing-inputs list.
