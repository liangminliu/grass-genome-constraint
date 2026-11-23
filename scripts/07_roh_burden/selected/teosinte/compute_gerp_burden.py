#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Teosinte genome-wide GERP burden: UNPOLARIZED ONLY.

Definition:
  raw ALT allele dosage from VCF genotype:
    0/0 -> 0
    0/1 -> 1
    1/1 -> 2

Input classes:
  - GERP4    : GERP > 4 SNPs
  - GERP_lt0 : GERP < 0 SNPs

Outputs:
  - Teosinte_GERP.site_AF.tsv
  - plot_A_site_spectrum.tsv
  - Teosinte_GERP.per_sample_burden.tsv
  - plot_BC_GERP4_burden.tsv
  - plot_DE_GERP4_over_lt0_ratio.tsv
  - plot_F_ROH_nonROH_burden.tsv
  - plot_GH_ROH_nonROH_ratio.tsv
  - Teosinte_GERP.diagnostics.tsv
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    sys.stderr.write("[ERROR] This script requires pysam. Install with: conda install -c bioconda pysam\n")
    raise

BASES = {"A", "C", "G", "T"}


def eprint(msg: str) -> None:
    sys.stderr.write(msg.rstrip() + "\n")


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def safe_div(num, den) -> float:
    try:
        den = float(den)
        num = float(num)
    except Exception:
        return float("nan")
    if den == 0 or math.isnan(den):
        return float("nan")
    return num / den


def fmt(x) -> str:
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.10g}"
    return str(x)


def normalize_chrom(chrom: str, mode: str = "auto") -> str:
    c = str(chrom)
    if mode == "exact":
        return c
    c = re.sub(r"^chr", "", c, flags=re.IGNORECASE)
    m = re.match(r"^0*([0-9]+)$", c)
    if m:
        return str(int(m.group(1)))
    return c


def parse_vcf_specs(specs: List[str]) -> List[Tuple[str, str, str]]:
    out = []
    for spec in specs:
        fields = spec.split(":", 2)
        if len(fields) != 3:
            raise ValueError(f"Bad --vcf spec: {spec}; expected class:source:path")
        cls, source, path = fields
        out.append((cls, source, path))
    return out


def read_froh(path: str) -> Dict[str, Dict[str, str]]:
    froh = {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get("sample") or row.get("IID") or row.get("id") or row.get("ID")
            if not sample:
                continue
            froh[sample] = row
            if row.get("IID"):
                froh[row["IID"]] = row
            if row.get("FID"):
                froh[row["FID"]] = row
    return froh


def find_col(header: List[str], candidates: List[str]) -> Optional[str]:
    lower = {h.lower(): h for h in header}
    for c in candidates:
        if c in header:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def read_roh_intervals(path: str, chrom_norm_mode: str = "auto") -> Dict[str, Dict[str, Tuple[List[int], List[int]]]]:
    """
    Load ROH intervals as:
      sample -> chrom -> (starts, ends)
    Coordinates are 1-based inclusive.
    """
    raw = defaultdict(lambda: defaultdict(list))

    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"ROH file has no header: {path}")

        header = reader.fieldnames
        sample_col = find_col(header, ["sample", "IID", "id", "ID"])
        chr_col = find_col(header, ["chr", "CHR", "chrom", "chromosome"])
        start_col = find_col(header, ["start_bp", "POS1", "BP1", "start", "START"])
        end_col = find_col(header, ["end_bp", "POS2", "BP2", "end", "END"])

        if not all([sample_col, chr_col, start_col, end_col]):
            raise RuntimeError(
                f"Could not detect ROH columns in {path}. "
                f"Need sample/IID, chr, start_bp/POS1, end_bp/POS2. Header={header}"
            )

        for row in reader:
            sample = str(row[sample_col])
            chrom = normalize_chrom(row[chr_col], chrom_norm_mode)
            try:
                st = int(float(row[start_col]))
                en = int(float(row[end_col]))
            except Exception:
                continue
            if en < st:
                st, en = en, st
            raw[sample][chrom].append((st, en))

    intervals = {}
    for sample, by_chr in raw.items():
        intervals[sample] = {}
        for chrom, vals in by_chr.items():
            vals.sort()
            merged = []
            for st, en in vals:
                if not merged or st > merged[-1][1] + 1:
                    merged.append((st, en))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], en))
            starts = [x[0] for x in merged]
            ends = [x[1] for x in merged]
            intervals[sample][chrom] = (starts, ends)

    return intervals


def pos_in_roh(intervals, sample: str, chrom: str, pos1: int) -> bool:
    by_chr = intervals.get(sample)
    if not by_chr:
        return False
    se = by_chr.get(chrom)
    if not se:
        return False
    starts, ends = se
    i = bisect.bisect_right(starts, pos1) - 1
    return i >= 0 and pos1 <= ends[i]


def init_counter():
    return {
        "called_sites": 0,
        "called_alleles": 0,
        "allele_count": 0,
        "hom_loci": 0,
        "het_loci": 0,
        "carrier_loci": 0,
        "additive_loci": 0.0,
    }


def update_counter(d: dict, ploidy: int, dosage: int) -> None:
    d["called_sites"] += 1
    d["called_alleles"] += ploidy
    d["allele_count"] += dosage

    if dosage == ploidy:
        d["hom_loci"] += 1

    if 0 < dosage < ploidy:
        d["het_loci"] += 1

    if dosage > 0:
        d["carrier_loci"] += 1

    # Diploid-equivalent additive load.
    # Diploid: hom=1, het=0.5.
    d["additive_loci"] += dosage / float(ploidy)


def extract_alt_dosage(gt):
    """
    Return (ploidy, ALT dosage) for biallelic GT.
    Missing or multi-allelic genotype is skipped.
    """
    if gt is None:
        return None
    alleles = list(gt)
    if len(alleles) == 0 or any(a is None for a in alleles):
        return None
    if any(a not in (0, 1) for a in alleles):
        return None
    ploidy = len(alleles)
    return ploidy, sum(1 for a in alleles if a == 1)


def make_af_bin(af: float) -> str:
    if math.isnan(af):
        return "NA"
    bins = [0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0000001]
    labels = [
        "0-0.01", "0.01-0.05", "0.05-0.10", "0.10-0.20",
        "0.20-0.30", "0.30-0.40", "0.40-0.50", "0.50-0.60",
        "0.60-0.70", "0.70-0.80", "0.80-0.90", "0.90-1.00"
    ]
    for i in range(len(bins) - 1):
        if bins[i] <= af < bins[i + 1]:
            return labels[i]
    if af == 1:
        return "0.90-1.00"
    return "NA"


def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--vcf", action="append", required=True, help="VCF spec class:source:path; repeatable")
    ap.add_argument("--froh", required=True, help="FROH.final.tsv")
    ap.add_argument("--roh", required=True, help="ROH.segments.classified.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chrom-normalize", choices=["auto", "exact"], default="auto")
    ap.add_argument("--skip-roh", action="store_true", help="Do not compute ROH/nonROH tables")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    vcf_specs = parse_vcf_specs(args.vcf)
    froh_map = read_froh(args.froh)
    eprint(f"[INFO] Loaded FROH records: {len(froh_map)} keys")

    roh_intervals = {} if args.skip_roh else read_roh_intervals(args.roh, args.chrom_normalize)
    eprint(f"[INFO] Loaded ROH interval samples: {len(roh_intervals)}")

    counts = defaultdict(init_counter)
    site_rows = []
    dropped = defaultdict(int)
    samples_seen = []
    sample_seen_set = set()

    def update_all(variant_class: str, chrom_norm: str, pos1: int, sample: str, ploidy: int, dosage: int) -> None:
        if sample not in sample_seen_set:
            sample_seen_set.add(sample)
            samples_seen.append(sample)

        update_counter(counts[(sample, variant_class, "ALL")], ploidy, dosage)

        if not args.skip_roh:
            region = "ROH" if pos_in_roh(roh_intervals, sample, chrom_norm, pos1) else "nonROH"
            update_counter(counts[(sample, variant_class, region)], ploidy, dosage)

    for variant_class, source, vcf_path in vcf_specs:
        eprint(f"[INFO] Processing {variant_class}:{source}: {vcf_path}")
        vf = pysam.VariantFile(vcf_path)

        for rec in vf.fetch():
            chrom = rec.chrom
            chrom_norm = normalize_chrom(chrom, args.chrom_normalize)
            pos1 = int(rec.pos)
            ref = rec.ref.upper() if rec.ref else ""
            alts = rec.alts or []

            if len(alts) != 1:
                dropped[(source, variant_class, "not_biallelic")] += 1
                continue

            alt = alts[0].upper()

            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
                dropped[(source, variant_class, "not_snp_acgt")] += 1
                continue

            acc = {
                "called_samples": 0,
                "called_alleles": 0,
                "allele_count": 0,
                "hom": 0,
                "het": 0,
                "carrier": 0,
            }

            for sample, call in rec.samples.items():
                x = extract_alt_dosage(call.get("GT"))
                if x is None:
                    continue

                ploidy, alt_dosage = x
                update_all(variant_class, chrom_norm, pos1, sample, ploidy, alt_dosage)

                acc["called_samples"] += 1
                acc["called_alleles"] += ploidy
                acc["allele_count"] += alt_dosage

                if alt_dosage == ploidy:
                    acc["hom"] += 1

                if 0 < alt_dosage < ploidy:
                    acc["het"] += 1

                if alt_dosage > 0:
                    acc["carrier"] += 1

            if acc["called_alleles"] > 0:
                af = safe_div(acc["allele_count"], acc["called_alleles"])
                site_rows.append({
                    "variant_class": variant_class,
                    "source": source,
                    "CHROM": chrom,
                    "POS": pos1,
                    "REF": ref,
                    "ALT": alt,
                    "target_allele": alt,
                    "target_type": "ALT",
                    "n_called_samples": acc["called_samples"],
                    "called_alleles": acc["called_alleles"],
                    "target_allele_count": acc["allele_count"],
                    "target_af": af,
                    "af_bin": make_af_bin(af),
                    "n_hom_loci": acc["hom"],
                    "n_het_loci": acc["het"],
                    "n_carrier_loci": acc["carrier"],
                })

        vf.close()

    # --------------------------------------------------------
    # 1. Site-level AF table
    # --------------------------------------------------------
    site_path = outdir / "Teosinte_GERP.site_AF.tsv"
    site_cols = [
        "variant_class", "source", "CHROM", "POS", "REF", "ALT",
        "target_allele", "target_type",
        "n_called_samples", "called_alleles", "target_allele_count",
        "target_af", "af_bin",
        "n_hom_loci", "n_het_loci", "n_carrier_loci",
    ]

    with open(site_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=site_cols)
        w.writeheader()
        for row in site_rows:
            w.writerow({k: fmt(row.get(k, "")) for k in site_cols})

    eprint(f"[INFO] Wrote {site_path}")

    # --------------------------------------------------------
    # 2. AF spectrum summary
    # --------------------------------------------------------
    spec_counts = defaultdict(int)
    spec_totals = defaultdict(int)

    for r in site_rows:
        key = (r["variant_class"], r["af_bin"])
        spec_counts[key] += 1
        spec_totals[r["variant_class"]] += 1

    spectrum_path = outdir / "plot_A_site_spectrum.tsv"
    with open(spectrum_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["variant_class", "af_bin", "n_sites", "pct_sites"])
        for (cls, afbin), n in sorted(spec_counts.items()):
            total = spec_totals[cls]
            w.writerow([cls, afbin, n, fmt(100 * safe_div(n, total))])

    eprint(f"[INFO] Wrote {spectrum_path}")

    # --------------------------------------------------------
    # 3. Per-sample burden table
    # --------------------------------------------------------
    burden_path = outdir / "Teosinte_GERP.per_sample_burden.tsv"
    burden_cols = [
        "sample", "variant_class", "region",
        "FROH", "n_roh", "sum_roh_bp", "Lauto_bp",
        "called_sites", "called_alleles",
        "target_allele_count", "hom_loci", "het_loci", "carrier_loci", "additive_loci",
        "allele_rate", "hom_loci_rate", "het_loci_rate", "carrier_loci_rate", "additive_load",
    ]

    with open(burden_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=burden_cols)
        w.writeheader()

        for (sample, cls, region), c in sorted(counts.items()):
            fr = froh_map.get(sample, {})
            row = {
                "sample": sample,
                "variant_class": cls,
                "region": region,
                "FROH": fr.get("FROH", "NA"),
                "n_roh": fr.get("n_roh", "NA"),
                "sum_roh_bp": fr.get("sum_roh_bp", "NA"),
                "Lauto_bp": fr.get("Lauto_bp", "NA"),
                "called_sites": c["called_sites"],
                "called_alleles": c["called_alleles"],
                "target_allele_count": c["allele_count"],
                "hom_loci": c["hom_loci"],
                "het_loci": c["het_loci"],
                "carrier_loci": c["carrier_loci"],
                "additive_loci": c["additive_loci"],
                "allele_rate": safe_div(c["allele_count"], c["called_alleles"]),
                "hom_loci_rate": safe_div(c["hom_loci"], c["called_sites"]),
                "het_loci_rate": safe_div(c["het_loci"], c["called_sites"]),
                "carrier_loci_rate": safe_div(c["carrier_loci"], c["called_sites"]),
                "additive_load": safe_div(c["additive_loci"], c["called_sites"]),
            }
            w.writerow({k: fmt(row.get(k, "")) for k in burden_cols})

    eprint(f"[INFO] Wrote {burden_path}")

    def get_count(sample: str, cls: str, region: str):
        return counts.get((sample, cls, region), init_counter())

    samples_seen = sorted(samples_seen)

    # --------------------------------------------------------
    # 4. B/C/D plot input: GERP > 4 burden
    # --------------------------------------------------------
    bc_path = outdir / "plot_BC_GERP4_burden.tsv"
    with open(bc_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "component", "FROH", "value_count", "value_rate", "called_sites"])

        for sample in samples_seen:
            c = get_count(sample, "GERP4", "ALL")
            fr = froh_map.get(sample, {}).get("FROH", "NA")

            comps = {
                "hom": (c["hom_loci"], safe_div(c["hom_loci"], c["called_sites"])),
                "het": (c["het_loci"], safe_div(c["het_loci"], c["called_sites"])),
                "carrier_all": (c["carrier_loci"], safe_div(c["carrier_loci"], c["called_sites"])),
                "additive": (c["additive_loci"], safe_div(c["additive_loci"], c["called_sites"])),
                "allele_count": (c["allele_count"], safe_div(c["allele_count"], c["called_alleles"])),
            }

            for comp, (val_count, val_rate) in comps.items():
                w.writerow([sample, comp, fr, fmt(val_count), fmt(val_rate), c["called_sites"]])

    eprint(f"[INFO] Wrote {bc_path}")

    # --------------------------------------------------------
    # 5. D/E ratio GERP > 4 / GERP < 0 for ALL region
    # --------------------------------------------------------
    states = [
        ("hom", "hom_loci", "hom_loci_rate"),
        ("het", "het_loci", "het_loci_rate"),
        ("carrier_all", "carrier_loci", "carrier_loci_rate"),
        ("additive", "additive_loci", "additive_load"),
        ("allele_count", "allele_count", "allele_rate"),
    ]

    ratio_path = outdir / "plot_DE_GERP4_over_lt0_ratio.tsv"
    with open(ratio_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sample", "state", "region", "FROH",
            "GERP4_value", "GERP_lt0_value", "ratio_count",
            "GERP4_rate", "GERP_lt0_rate", "ratio_rate",
        ])

        for sample in samples_seen:
            fr = froh_map.get(sample, {}).get("FROH", "NA")
            c4 = get_count(sample, "GERP4", "ALL")
            c0 = get_count(sample, "GERP_lt0", "ALL")

            metrics4 = {
                "hom_loci_rate": safe_div(c4["hom_loci"], c4["called_sites"]),
                "het_loci_rate": safe_div(c4["het_loci"], c4["called_sites"]),
                "carrier_loci_rate": safe_div(c4["carrier_loci"], c4["called_sites"]),
                "additive_load": safe_div(c4["additive_loci"], c4["called_sites"]),
                "allele_rate": safe_div(c4["allele_count"], c4["called_alleles"]),
            }

            metrics0 = {
                "hom_loci_rate": safe_div(c0["hom_loci"], c0["called_sites"]),
                "het_loci_rate": safe_div(c0["het_loci"], c0["called_sites"]),
                "carrier_loci_rate": safe_div(c0["carrier_loci"], c0["called_sites"]),
                "additive_load": safe_div(c0["additive_loci"], c0["called_sites"]),
                "allele_rate": safe_div(c0["allele_count"], c0["called_alleles"]),
            }

            for state, count_key, rate_key in states:
                v4 = c4["allele_count"] if count_key == "allele_count" else c4[count_key]
                v0 = c0["allele_count"] if count_key == "allele_count" else c0[count_key]
                r4 = metrics4[rate_key]
                r0 = metrics0[rate_key]

                w.writerow([
                    sample, state, "ALL", fr,
                    fmt(v4), fmt(v0), fmt(safe_div(v4, v0)),
                    fmt(r4), fmt(r0), fmt(safe_div(r4, r0)),
                ])

    eprint(f"[INFO] Wrote {ratio_path}")

    # --------------------------------------------------------
    # 6. ROH vs nonROH burden and ratio
    # --------------------------------------------------------
    roh_burden_path = outdir / "plot_F_ROH_nonROH_burden.tsv"
    roh_ratio_path = outdir / "plot_GH_ROH_nonROH_ratio.tsv"

    with open(roh_burden_path, "w", newline="") as fb, open(roh_ratio_path, "w", newline="") as frh:
        wb = csv.writer(fb, delimiter="\t")
        wr = csv.writer(frh, delimiter="\t")

        wb.writerow([
            "sample", "variant_class", "region", "component", "FROH",
            "value_count", "value_rate", "called_sites",
        ])

        wr.writerow([
            "sample", "state", "region", "FROH",
            "GERP4_value", "GERP_lt0_value", "ratio_count",
            "GERP4_rate", "GERP_lt0_rate", "ratio_rate",
        ])

        for sample in samples_seen:
            fr = froh_map.get(sample, {}).get("FROH", "NA")

            for region in ["ROH", "nonROH"]:
                for cls in ["GERP4", "GERP_lt0"]:
                    c = get_count(sample, cls, region)

                    comps = {
                        "hom": (c["hom_loci"], safe_div(c["hom_loci"], c["called_sites"])),
                        "het": (c["het_loci"], safe_div(c["het_loci"], c["called_sites"])),
                        "carrier_all": (c["carrier_loci"], safe_div(c["carrier_loci"], c["called_sites"])),
                        "additive": (c["additive_loci"], safe_div(c["additive_loci"], c["called_sites"])),
                        "allele_count": (c["allele_count"], safe_div(c["allele_count"], c["called_alleles"])),
                    }

                    for comp, (val_count, val_rate) in comps.items():
                        wb.writerow([sample, cls, region, comp, fr, fmt(val_count), fmt(val_rate), c["called_sites"]])

                c4 = get_count(sample, "GERP4", region)
                c0 = get_count(sample, "GERP_lt0", region)

                rate_lookup4 = {
                    "hom_loci_rate": safe_div(c4["hom_loci"], c4["called_sites"]),
                    "het_loci_rate": safe_div(c4["het_loci"], c4["called_sites"]),
                    "carrier_loci_rate": safe_div(c4["carrier_loci"], c4["called_sites"]),
                    "additive_load": safe_div(c4["additive_loci"], c4["called_sites"]),
                    "allele_rate": safe_div(c4["allele_count"], c4["called_alleles"]),
                }

                rate_lookup0 = {
                    "hom_loci_rate": safe_div(c0["hom_loci"], c0["called_sites"]),
                    "het_loci_rate": safe_div(c0["het_loci"], c0["called_sites"]),
                    "carrier_loci_rate": safe_div(c0["carrier_loci"], c0["called_sites"]),
                    "additive_load": safe_div(c0["additive_loci"], c0["called_sites"]),
                    "allele_rate": safe_div(c0["allele_count"], c0["called_alleles"]),
                }

                for state, count_key, rate_key in states:
                    v4 = c4["allele_count"] if count_key == "allele_count" else c4[count_key]
                    v0 = c0["allele_count"] if count_key == "allele_count" else c0[count_key]
                    r4 = rate_lookup4[rate_key]
                    r0 = rate_lookup0[rate_key]

                    wr.writerow([
                        sample, state, region, fr,
                        fmt(v4), fmt(v0), fmt(safe_div(v4, v0)),
                        fmt(r4), fmt(r0), fmt(safe_div(r4, r0)),
                    ])

    eprint(f"[INFO] Wrote {roh_burden_path}")
    eprint(f"[INFO] Wrote {roh_ratio_path}")

    # --------------------------------------------------------
    # 7. Diagnostics
    # --------------------------------------------------------
    diag_path = outdir / "Teosinte_GERP.diagnostics.tsv"
    with open(diag_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source", "variant_class", "reason", "n"])
        for (source, cls, reason), n in sorted(dropped.items()):
            w.writerow([source, cls, reason, n])

    eprint(f"[INFO] Wrote {diag_path}")
    eprint("[DONE] Teosinte unpolarized calculation finished.")


if __name__ == "__main__":
    main()
