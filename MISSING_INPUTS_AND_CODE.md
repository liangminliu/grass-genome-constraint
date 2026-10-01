# Remaining Methods inputs and code to provide or confirm

| Methods step | Needed to complete the released workflow |
| --- | --- |
| Grass phylogeny | 97-genome FASTA/accession manifest, production ROADIES log and final tree hash; final treePL calibration/configuration and dated-tree output if treePL code is to be released. |
| Genome-wide synteny | Final GeneTribe 3 × 93 pair manifest and launch script, reference/query assembly versions, block-retention rule and final collinear-block table; JCVI command and figure input for the chromosome view. |
| Whole-genome alignment | Exact existing ROADIES tree, focal YAML/Snakefile revision, 93-query lists for each reference, three `roast.maf` hashes and the BED/coverage tables used for Fig. 2. The release MSA driver depends on the original pipeline checkout, which is not bundled. |
| GERP constraint | Final fourfold-site alignment, neutral tree/phyloFit model, masks, GERP++ command and per-base score files for each reference; category BED definitions and curated tables used by the pie/histogram scripts. |
| Population SNP filtering | Bamboo and teosinte final mapping→GATK→KING→filter commands and sample sheets; final VCF hashes. The supplied A. tauschii list specifies 141 IDs; provide any separate 140-ID non-ROH analysis input if used. |
| Deleterious annotation | SIFT 4G and SnpEff database builds/output hashes, explicit SIFT <0.05 extraction command, GERP/SIFT/functional-effect join keys, allele-frequency bins and denominator definitions. |
| ROH/FROH and burden | Species-specific `Lauto.from_filtered_SNPs.tsv` inputs, final PLINK `.hom` and FROH outputs, exact burden/heterozygosity joins and the final Fig. 6/S8 table-to-panel map. Confirm whether `--homozyg-het 1` was explicit or the PLINK default in the production log. |
| Synteny-based gene comparisons | Exact five SynPan `.prot` and `.bed` inputs, final `grass.list.SG` hash, gene-class assignment script/output, 5,108 matched-group source table, and each Fig. 7/S9/S10 test-to-panel link. |

These records can be added one step at a time. The current public tree intentionally has no `source_candidates/` or `drafts/` directories.
