# Heterozygosity, ROH and burden

The species-specific `run_roh_froh.sh` scripts call ROH/FROH with gap100/kb100 and five missing sites per window. `selected/common/run_gerp_burden.sh` runs the three recovered GERP burden calculators using the final ROH and FROH inputs supplied at run time. It produces per-sample tables for S6 and S8; the associated plotting programs are in this directory.

The recovered original burden launchers referenced earlier `gap1000` ROH paths, so they are not used by the shared launcher. Confirm the ROH/FROH source underlying the published S3/S6 tables before claiming an exact end-to-end rerun.
