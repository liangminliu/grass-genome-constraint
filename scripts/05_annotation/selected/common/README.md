# Shared SIFT 4G and SnpEff commands

Run from the species analysis working directory. Pass absolute paths to the workflow and settings file:

```bash
REPO=/path/to/grass-genome-constraint
FLOW="$REPO/scripts/05_annotation/selected/common/workflow.sh"
SETTINGS="$REPO/scripts/05_annotation/selected/atau/settings.sh"  # or ped/zmay
bash "$FLOW" "$SETTINGS" sift-build
# Wait for the chromosome array, then:
bash "$FLOW" "$SETTINGS" sift-merge
bash "$FLOW" "$SETTINGS" sift-annotate
bash "$FLOW" "$SETTINGS" snpeff-build
# Wait for the database job, then:
bash "$FLOW" "$SETTINGS" snpeff-all
bash "$FLOW" "$SETTINGS" snpeff-gerp
```

`snpeff-sift` annotates the existing SIFT-derived VCF specified for Zmay or Atau. Bamboo's separate C/D SIFT classes await the exact final class VCFs and are listed in the repository's missing-inputs file. The workflow submits heavy jobs with SLURM and runs the symlink merge locally; wait for dependent jobs before the next action. `settings.sh` contains source-cluster paths that must exist or be updated when installing elsewhere. The source jobs used `-ud 5000` for SnpEff. The common merge checks basename collisions and regenerates `all_prot.fasta` on each run instead of appending duplicate sequences.
