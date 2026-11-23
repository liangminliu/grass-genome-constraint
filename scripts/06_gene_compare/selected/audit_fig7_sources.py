#!/usr/bin/env python3
"""Validate Fig. 7/S9/S10 gene-level inputs and write panel-to-test mapping."""
import argparse
import csv
import json
from pathlib import Path


MAP = [
    ('Fig7A;S9A', 'strict', 'tables/Strict_core_SG_GERP_CDS_wide_all_thresholds.tsv', 'stats/Strict_core_GERP_CDS_paired_stats.tsv', 'cds_conserved_pct', 'paired SG', 'gt2,gt4,gt6'),
    ('Fig7B;S10A', 'class', 'tables/GERP_CDS_constrained_pct_all_species_thresholds.tsv', 'stats/GERP_CDS_constrained_pct_stats_all_species_thresholds.tsv', 'cds_conserved_pct', 'individual gene', 'gt2,gt4,gt6'),
    ('Fig7C;S10B', 'class', 'tables/SNP_metrics_gt4_ALL.tsv', 'stats/SNP_stats_gt4_ALL.tsv', 'cds_snp_density', 'individual gene', 'gt4 input'),
    ('Fig7D;S10C', 'class', 'tables/SNP_metrics_all_species_thresholds_data_types.tsv', 'stats/SNP_stats_all_species_thresholds_data_types.tsv', 'cds_deleterious_density', 'individual gene', 'gt2,gt4,gt6'),
    ('Fig7E;S10D', 'class', 'tables/SNP_metrics_all_species_thresholds_data_types.tsv', 'stats/SNP_stats_all_species_thresholds_data_types.tsv', 'cds_deleterious_ratio', 'SNP-positive individual gene', 'gt2,gt4,gt6'),
    ('Fig7F;S9B', 'strict', 'tables/Strict_core_SG_SNP_deleterious_wide_all_thresholds.tsv', 'stats/Strict_core_SNP_deleterious_stats.tsv', 'cds_snp_density', 'paired SG', 'gt4 input'),
    ('Fig7G;S9D', 'strict', 'tables/Strict_core_SG_SNP_deleterious_wide_all_thresholds.tsv', 'stats/Strict_core_SNP_deleterious_stats.tsv', 'cds_deleterious_density', 'paired SG', 'gt2,gt4,gt6'),
    ('Fig7H;S9E', 'strict', 'tables/Strict_core_SG_SNP_deleterious_wide_all_thresholds.tsv', 'stats/Strict_core_SNP_deleterious_stats.tsv', 'cds_deleterious_ratio', 'SNP-positive individual gene (ALL); paired SG (matched subset)', 'gt2,gt4,gt6'),
    ('S9C', 'strict', 'tables/Strict_core_SNP0_SG_proportion.tsv', 'stats/Strict_core_SNP0_CochranQ_exact_McNemar.tsv', 'SNP0 proportion', 'paired binary SG', 'gt2,gt4,gt6'),
]


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--strict-dir', type=Path, required=True)
    p.add_argument('--class-dir', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    a = p.parse_args()
    roots = {'strict': a.strict_dir, 'class': a.class_dir}
    absent = sorted({str(roots[k] / rel) for _, k, t, s, *_ in MAP for rel in (t, s) if not (roots[k] / rel).is_file()})
    if absent:
        raise SystemExit('Missing source files:\n' + '\n'.join(absent))
    wide = a.strict_dir / 'tables/Strict_core_SG_GERP_CDS_wide_all_thresholds.tsv'
    group_ids = {r['SG'] for r in rows(wide) if r['threshold'] == 'gt4'}
    if len(group_ids) != 5108:
        raise SystemExit(f'Expected 5,108 strict-core groups at gt4; found {len(group_ids)}')
    subset = list(rows(a.strict_dir / 'tables/Strict_core_SNPpos_subset_definition.tsv'))
    if len(subset) != 1 or int(subset[0]['n_canonical_SNP_positive_strict_SG']) != 668:
        raise SystemExit('Expected 668 matched SNP-positive strict-core groups')
    snp_wide = a.strict_dir / 'tables/Strict_core_SG_SNP_deleterious_wide_all_thresholds.tsv'
    zero_filled = {r['SG'] for r in rows(snp_wide)
                   if r['threshold'] == 'gt4' and r.get('Atau_absent_SNP_row_zero_filled') == 'TRUE'}
    if zero_filled != {'SG0039462', 'SG0039466'}:
        raise SystemExit(f'Unexpected A. tauschii zero-filled groups: {sorted(zero_filled)}')
    paired = list(rows(a.strict_dir / 'stats/Strict_core_GERP_CDS_paired_stats.tsv'))
    if any(int(r['n_pairs']) != 5108 for r in paired):
        raise SystemExit('Paired GERP test does not use 5,108 groups')
    maize_atau = [r for r in paired if r['threshold'] == 'gt2' and r['Comparison'] == 'Z.mays - A.tauschii']
    if len(maize_atau) != 1 or abs(float(maize_atau[0]['P.adj']) - 0.0972031304835916) > 1e-10:
        raise SystemExit('GERP >2 maize–A. tauschii bracket source differs from audited value')
    a.outdir.mkdir(parents=True, exist_ok=True)
    mapping = a.outdir / 'panel_to_source_test.tsv'
    with mapping.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('panel', 'source_set', 'gene_table', 'test_table', 'metric', 'test_unit', 'thresholds'))
        writer.writerows(MAP)
    audit = {'strict_core_groups_gt4': len(group_ids), 'matched_snp_positive_groups': 668,
             'aegilops_zero_filled_groups': sorted(zero_filled),
             'gerp_gt2_maize_atau_bh_p': float(maize_atau[0]['P.adj']),
             'source_files_checked': len({(k, rel) for _, k, t, s, *_ in MAP for rel in (t, s)})}
    (a.outdir / 'source_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
