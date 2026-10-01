#!/usr/bin/env bash
#BSUB -J plot_gerp_pubstyle
#BSUB -q Q64C1T_X4
#BSUB -n 10
#BSUB -o logs/plot_gerp_pubstyle.%J.out
#BSUB -e logs/plot_gerp_pubstyle.%J.err

# Run from a directory containing the three input TSV files.
# Create logs/ before submitting this job so LSF can open its log paths.
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
plot_script="$script_dir/plot_gerp_mbe_final.py"
for input in ALL_gerp.tsv C_CDS_gerp.tsv D_CDS_gerp.tsv; do
  [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
[[ -s "$plot_script" ]] || { echo "Missing plot script: $plot_script" >&2; exit 1; }
mkdir -p logs

echo "Starting GERP histogram: $(date)"
echo "Host: $(hostname); working directory: $(pwd)"
python3 "$plot_script" \
  --all_file ALL_gerp.tsv \
  --c_cds_file C_CDS_gerp.tsv \
  --d_cds_file D_CDS_gerp.tsv \
  --bin_size 0.4 \
  --output_prefix ped_gerp_mbe_final_0.4v5 \
  --annotate_overflow \
  --show_break_note \
  --fig_width 4.8 \
  --fig_height 4.2
echo "Finished GERP histogram: $(date)"
