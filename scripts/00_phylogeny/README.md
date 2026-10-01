# Grass phylogeny

Run `selected/run_roadies_97sp_accurate_v2.sh` with `config/config_97sp_accurate.yaml` for 97 genome inputs (94 grasses plus three outgroups). The supplied job uses ROADIES accurate mode, continuity 85, support threshold 0.95 and `GENE_COUNT: 16000` initial sampled fragments. Update the cluster paths before submission. The job clears its configured output directories at startup; save any existing output first.

The final treePL invocation and calibration configuration remain pending, so this release stops at the confirmed ROADIES job. See the root missing-inputs list.
