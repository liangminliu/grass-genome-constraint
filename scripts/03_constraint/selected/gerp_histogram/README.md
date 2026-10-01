# Bamboo genome-wide and CDS GERP histogram

`plot_gerp_mbe_final.py` is the exact author-supplied Python plotting code. It reads the first column of `ALL_gerp.tsv`, `C_CDS_gerp.tsv`, and `D_CDS_gerp.tsv` in chunks; non-numeric rows are ignored. It bins GERP scores at width 0.4, combines the two bamboo CDS components in a stacked display, uses a broken y-axis, and writes PDF and 300-dpi PNG outputs.

`run_gerp_histogram.lsf.sh` is the corrected LSF wrapper. It checks the Python script actually invoked, preserves the supplied plotting parameters, and exits on failure. The author assigned this distribution to Fig. 4A; the original wrapper's internal panel label `B` is suppressed in this release so the final multi-panel assembly can apply its own label. From the input directory, create `logs/` before submitting the job. Outputs are `ped_gerp_mbe_final_0.4v5.pdf` and `.png`.

**Figure mapping:** The current manuscript captions Fig. 3A as a chromosome-track view and Fig. 4A as the three-species GERP-score distribution. This code draws only the bamboo histogram and labels it `B`; it is a component candidate for the Fig. 4A distribution, not a stand-alone reproduction of Fig. 3A. Confirm the final assembled figure and species-panel labels against the manuscript image before changing this job's label or citing it as the final panel.

The script needs Python, NumPy, pandas and Matplotlib. The large GERP TSV inputs are excluded from this release; their final hashes and the equivalent maize/goatgrass histogram jobs remain to be linked to the assembled figure.
