#!/usr/bin/env python3
"""Build the five 500-kb tracks used by Fig. 3 and Fig. S1.

All interval inputs are BED-style, zero-based and half-open. GERP and pi files
are four-column bedGraphs (chrom, start, end, value). Gene and TE files contain
at least the first three BED columns. The SNP file is a biallelic VCF/VCF.GZ.
"""
import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.strip() and not line.startswith('#'):
                yield line.rstrip('\n').split('\t')


def sizes(path):
    result = {}
    for row in lines(path):
        if len(row) < 2:
            continue
        chrom, length = row[0], int(row[1])
        if chrom in result or length <= 0:
            raise ValueError(f'Duplicate/invalid chromosome: {chrom}')
        result[chrom] = length
    if not result:
        raise ValueError('Chromosome-size table is empty')
    return result


def overlap_windows(chrom, start, end, lengths, width):
    if chrom not in lengths:
        return
    if start < 0 or end <= start or end > lengths[chrom]:
        raise ValueError(f'Invalid interval: {chrom}:{start}-{end}')
    for i in range(start // width, (end - 1) // width + 1):
        left, right = i * width, min((i + 1) * width, lengths[chrom])
        yield i, min(end, right) - max(start, left)


def add_bedgraph(path, lengths, width, numerators, denominators, label):
    previous = {}
    for row in lines(path):
        if len(row) < 4:
            raise ValueError(f'{label}: expected four bedGraph columns')
        chrom, start, end = row[0], int(row[1]), int(row[2])
        value = float(row[3])
        if not math.isfinite(value):
            continue
        if chrom not in lengths:
            continue
        if start < previous.get(chrom, 0):
            raise ValueError(f'{label} bedGraph overlaps or is unsorted on {chrom}')
        previous[chrom] = end
        for i, n in overlap_windows(chrom, start, end, lengths, width):
            numerators[(chrom, i)] += value * n
            denominators[(chrom, i)] += n


def add_merged_bed(path, lengths, width, covered):
    spans = defaultdict(list)
    for row in lines(path):
        if len(row) < 3:
            raise ValueError('BED requires at least three columns')
        chrom, start, end = row[0], int(row[1]), int(row[2])
        if chrom in lengths:
            if start < 0 or end <= start or end > lengths[chrom]:
                raise ValueError(f'Invalid BED interval: {chrom}:{start}-{end}')
            spans[chrom].append((start, end))
    for chrom, intervals in spans.items():
        intervals.sort()
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        for start, end in merged:
            for i, n in overlap_windows(chrom, start, end, lengths, width):
                covered[(chrom, i)] += n


def add_vcf(path, lengths, width, counts):
    seen = set()
    for row in lines(path):
        if len(row) < 5 or row[0] not in lengths:
            continue
        chrom, pos, ref, alt = row[0], int(row[1]), row[3], row[4]
        if len(ref) != 1 or len(alt) != 1 or alt not in 'ACGT' or ref not in 'ACGT':
            continue
        if not 1 <= pos <= lengths[chrom]:
            raise ValueError(f'VCF position outside chromosome: {chrom}:{pos}')
        key = (chrom, pos, ref, alt)
        if key not in seen:
            counts[(chrom, (pos - 1) // width)] += 1
            seen.add(key)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('chrom-sizes', 'gerp-bedgraph', 'genes-bed', 'te-bed', 'snp-vcf', 'pi-bedgraph', 'output'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--window-bp', type=int, default=500000)
    a = p.parse_args()
    if a.window_bp != 500000:
        raise ValueError('The manuscript Fig. 3/S1 window size is 500000 bp')
    lengths = sizes(a.chrom_sizes)
    gerp_sum, gerp_n = defaultdict(float), defaultdict(int)
    pi_sum, pi_n = defaultdict(float), defaultdict(int)
    genes, te, snps = defaultdict(int), defaultdict(int), defaultdict(int)
    add_bedgraph(a.gerp_bedgraph, lengths, a.window_bp, gerp_sum, gerp_n, 'GERP')
    add_bedgraph(a.pi_bedgraph, lengths, a.window_bp, pi_sum, pi_n, 'pi')
    add_merged_bed(a.genes_bed, lengths, a.window_bp, genes)
    add_merged_bed(a.te_bed, lengths, a.window_bp, te)
    add_vcf(a.snp_vcf, lengths, a.window_bp, snps)
    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('chrom', 'start', 'end', 'gerp_mean', 'gerp_scored_bp',
                         'gene_fraction', 'te_fraction', 'snp_count', 'snp_per_kb',
                         'pi_mean', 'pi_scored_bp'))
        for chrom, length in lengths.items():
            for i in range((length + a.window_bp - 1) // a.window_bp):
                start, end = i * a.window_bp, min((i + 1) * a.window_bp, length)
                key = (chrom, i)
                span = end - start
                writer.writerow((chrom, start, end,
                    f'{gerp_sum[key] / gerp_n[key]:.8g}' if gerp_n[key] else 'NA',
                    gerp_n[key], f'{genes[key] / span:.8g}', f'{te[key] / span:.8g}',
                    snps[key], f'{snps[key] * 1000 / span:.8g}',
                    f'{pi_sum[key] / pi_n[key]:.8g}' if pi_n[key] else 'NA', pi_n[key]))


if __name__ == '__main__':
    main()
