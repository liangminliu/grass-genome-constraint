#!/usr/bin/env python3
"""Audit final S6 sheets against production GERP>4 burden tables for Fig. S8."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PANELS = (
    ('bamboo', 'S6C_bamboo', 63),
    ('teosinte', 'S6B_teosinte', 90),
    ('Aegilops', 'S6A_Aegilops', 140),
)
METRICS = ('FROH', 'het_loci', 'hom_loci', 'target_allele_count')


def select(df):
    df = df.copy()
    klass = df['variant_class'].astype(str).str.replace(' ', '', regex=False)
    keep = df['region'].eq('ALL') & klass.isin(('GERP4', 'GERP>4'))
    if 'mode' in df:
        keep &= df['mode'].eq('unpolarized')
    df = df.loc[keep].rename(columns={'alt_allele_count': 'target_allele_count'})
    if df['sample'].duplicated().any():
        raise ValueError('Duplicate sample after ALL/GERP>4 selection')
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--s6-workbook', type=Path, required=True)
    for name, _, _ in PANELS:
        p.add_argument('--' + name.lower() + '-burden', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    assembled, audit = [], []
    for name, sheet, expected in PANELS:
        source = getattr(args, name.lower() + '_burden')
        approved = select(pd.read_excel(args.s6_workbook, sheet_name=sheet, header=1))
        raw = select(pd.read_csv(source, sep='\t'))
        if len(approved) != expected:
            raise ValueError(f'{sheet}: expected {expected} samples, found {len(approved)}')
        if not set(approved['sample']).issubset(set(raw['sample'])):
            raise ValueError(f'{name}: S6 sample missing from raw burden table')
        joined = approved.merge(raw, on='sample', suffixes=('_s6', '_raw'), validate='one_to_one')
        if len(joined) != expected:
            raise ValueError(f'{name}: failed one-to-one source match')
        deltas = {}
        for metric in METRICS:
            left = pd.to_numeric(joined[metric + '_s6'], errors='raise')
            right = pd.to_numeric(joined[metric + '_raw'], errors='raise')
            delta = float(np.max(np.abs(left - right)))
            if not np.isfinite(delta) or delta > 1e-7:
                raise ValueError(f'{name}: {metric} differs from raw data by {delta}')
            deltas[metric] = delta
        if not np.allclose(approved['target_allele_count'],
                           2 * approved['hom_loci'] + approved['het_loci'], atol=0, rtol=0):
            raise ValueError(f'{name}: allele dosage does not match genotype-state counts')
        part = approved[['sample', *METRICS]].copy()
        part.insert(0, 'species', name)
        assembled.append(part)
        audit.append({'species': name, 'n_s6': expected, 'n_raw': len(raw),
                      'excluded_from_s6': sorted(set(raw['sample']) - set(approved['sample'])),
                      'max_absolute_difference': deltas,
                      'raw_sha256': hashlib.sha256(source.read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(assembled, ignore_index=True).to_csv(args.output, sep='\t', index=False)
    audit_path = args.output.with_suffix('.audit.json')
    audit_path.write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(f'Validated {sum(x[2] for x in PANELS)} individuals; wrote {args.output} and {audit_path}')


if __name__ == '__main__':
    main()
