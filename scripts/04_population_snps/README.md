# Population SNP filtering

- `selected/common/`: reads to gVCF, joint SNP calling, relatedness filtering and final population VCF filtering for bamboo and maize/teosinte. These upstream commands were reconstructed from Methods.
- `selected/atau/`: supplied *A. tauschii* VCF filtering and sample selection scripts.

Provide each reference genome and sequencing inputs. The *A. tauschii* retained-sample list is in `selected/atau/retained_samples.txt`; pass its path through `SAMPLES` when running the extraction script.

