# Population SNP filtering

`selected/common/` contains the bamboo and teosinte SNP workflow. The retained-sample lists in `selected/ped/` and `selected/zmay/` were exported from final Supplementary Table S3 and match the sample sets in Table S6. Pass the corresponding list to `04_final_population_filter.sh` at run time. `extract_retained_lists_from_s3.py` regenerates both lists from the workbook.

`selected/atau/` contains the supplied *A. tauschii* filtering script and retained-sample list. Pass the list path through `SAMPLES` when running its extraction script.
