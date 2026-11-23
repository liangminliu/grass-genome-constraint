#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
import re
from bisect import bisect_right
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_sample(x):
    return re.sub(r"\s+", "", str(x).strip())


def norm_chr(x):
    x = str(x).strip()
    x = re.sub(r"^chr", "", x, flags=re.I)
    x = x.upper()
    x = re.sub(r"\.0$", "", x)
    if re.fullmatch(r"[1-7]", x):
        x = x + "D"
    return x


def parse_exclude(x):
    if not x or x.lower() in {"none", "na", "null"}:
        return set()
    return {norm_sample(i) for i in re.split(r"[,; \t]+", x) if i.strip()}


def parse_gt(field):
    gt = field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    if sep is None:
        return None
    arr = gt.split(sep)
    if any(a == "." for a in arr):
        return None
    try:
        return [int(a) for a in arr]
    except Exception:
        return None


def is_biallelic_snp(ref, alt):
    return (
        len(ref) == 1 and len(alt) == 1 and "," not in alt and
        ref.upper() in {"A", "C", "G", "T"} and
        alt.upper() in {"A", "C", "G", "T"}
    )


def pick_col(cols, candidates, required=True):
    for c in candidates:
        if c in cols:
            return c
    if required:
        raise ValueError(f"Cannot find column. Tried: {candidates}; current columns: {list(cols)}")
    return None


def read_lauto_total(path):
    total = None
    with open(path, "rt", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r.get("chr") == "TOTAL":
                total = int(float(r["Lauto_bp"]))
    if total is None:
        raise ValueError(f"Cannot find TOTAL Lauto_bp in {path}")
    return total


def read_froh(path):
    out = {}
    with open(path, "rt", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = list(rd)

    if not rows:
        return out

    sid_col = pick_col(rows[0].keys(), ["sample", "IID", "Sample", "ID"], required=True)

    for r in rows:
        sid = norm_sample(r.get(sid_col, ""))
        if sid:
            out[sid] = r

    return out


def read_roh_segments(path, exclude):
    raw = defaultdict(lambda: defaultdict(list))

    with open(path, "rt", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames

        iid_c = pick_col(cols, ["sample", "IID", "Sample", "ID"])
        chr_c = pick_col(cols, ["chr", "CHR", "chrom", "chromosome"])
        s_c = pick_col(cols, ["start_bp", "POS1", "start"])
        e_c = pick_col(cols, ["end_bp", "POS2", "end"])

        for r in rd:
            iid = norm_sample(r[iid_c])
            if iid in exclude:
                continue

            chrom = norm_chr(r[chr_c])
            s = int(float(r[s_c]))
            e = int(float(r[e_c]))

            if e < s:
                s, e = e, s

            raw[iid][chrom].append((s, e))

    merged = {}

    for sid, chrom_map in raw.items():
        merged[sid] = {}
        for chrom, arr in chrom_map.items():
            arr = sorted(arr)
            if not arr:
                continue

            out = []
            cs, ce = arr[0]

            for s, e in arr[1:]:
                if s <= ce + 1:
                    ce = max(ce, e)
                else:
                    out.append((cs, ce))
                    cs, ce = s, e

            out.append((cs, ce))
            merged[sid][chrom] = out

    return merged


def build_interval_index(rohs):
    idx = {}
    for sid, chrom_map in rohs.items():
        idx[sid] = {}
        for chrom, arr in chrom_map.items():
            starts = [s for s, e in arr]
            ends = [e for s, e in arr]
            idx[sid][chrom] = (starts, ends)
    return idx


def in_roh(interval_index, sample, chrom, pos):
    chrom = norm_chr(chrom)
    d = interval_index.get(sample, {})
    if chrom not in d:
        return False
    starts, ends = d[chrom]
    i = bisect_right(starts, pos) - 1
    return i >= 0 and pos <= ends[i]


def init_count():
    return {
        "called_sites": 0,
        "called_alleles": 0,
        "alt_allele_count": 0,
        "hom_loci": 0,
        "het_loci": 0,
    }


def update_count(d, alleles, ploidy):
    d["called_sites"] += 1
    d["called_alleles"] += len(alleles)

    alt_n = sum(1 for a in alleles if a == 1)

    d["alt_allele_count"] += alt_n

    if alt_n == ploidy:
        d["hom_loci"] += 1
    elif 0 < alt_n < ploidy:
        d["het_loci"] += 1


def safe_div(a, b):
    return a / b if b else 0.0


def parse_vcf(vcf, class_label, exclude, ploidy, interval_index, counts, site_rows):
    with open_text(vcf) as fh:
        samples = []
        keep_i = []

        for line in fh:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                h = line.rstrip("\n").split("\t")
                samples_all = [norm_sample(x) for x in h[9:]]

                for i, s in enumerate(samples_all):
                    if s not in exclude:
                        samples.append(s)
                        keep_i.append(i)
                break

        if not samples:
            raise ValueError(f"No samples found in VCF: {vcf}")

        for s in samples:
            for region in ["ALL", "ROH", "Non-ROH"]:
                counts[(s, class_label, region)]

        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            p = line.rstrip("\n").split("\t")
            chrom, pos_s, vid, ref, alt = p[:5]

            if not is_biallelic_snp(ref, alt):
                continue

            chrom_norm = norm_chr(chrom)
            pos = int(pos_s)

            ac = 0
            an = 0
            n_called = 0
            hom = 0
            het = 0

            sample_fields = p[9:]

            for sample, orig_i in zip(samples, keep_i):
                if orig_i >= len(sample_fields):
                    continue

                gt = parse_gt(sample_fields[orig_i])

                if gt is None:
                    continue

                if len(gt) != ploidy:
                    continue

                if any(a not in {0, 1} for a in gt):
                    continue

                n_called += 1
                alt_n = sum(1 for a in gt if a == 1)
                ac += alt_n
                an += len(gt)

                if alt_n == ploidy:
                    hom += 1
                elif 0 < alt_n < ploidy:
                    het += 1

                inside = in_roh(interval_index, sample, chrom_norm, pos)

                update_count(counts[(sample, class_label, "ALL")], gt, ploidy)

                if inside:
                    update_count(counts[(sample, class_label, "ROH")], gt, ploidy)
                else:
                    update_count(counts[(sample, class_label, "Non-ROH")], gt, ploidy)

            if an > 0:
                site_rows.append({
                    "variant_class": class_label,
                    "chr": chrom_norm,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "AC": ac,
                    "AN": an,
                    "AF": safe_div(ac, an),
                    "n_called_samples": n_called,
                    "hom_loci": hom,
                    "het_loci": het,
                })


def write_outputs(out_prefix, counts, site_rows, froh, lgenome):
    site_file = out_prefix + ".site_AF.tsv"
    burden_file = out_prefix + ".per_sample_burden.tsv"
    ratio_file = out_prefix + ".ratio_metrics.tsv"

    with open(site_file, "wt", newline="") as fh:
        fieldnames = [
            "variant_class", "chr", "pos", "ref", "alt",
            "AC", "AN", "AF", "n_called_samples", "hom_loci", "het_loci"
        ]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in site_rows:
            w.writerow(r)

    burden_rows = []
    ratio_rows = []

    for (sample, cls, region), c in sorted(counts.items()):
        frow = froh.get(sample, {})

        FROH = frow.get("FROH", frow.get("Froh", frow.get("froh", "")))
        SROH = frow.get("sum_roh_bp", frow.get("SROH_bp", ""))
        NROH = frow.get("n_roh", frow.get("NROH", ""))

        row = {
            "sample": sample,
            "variant_class": cls,
            "region": region,
            "Lgenome_bp": lgenome,
            "FROH": FROH,
            "SROH_bp": SROH,
            "NROH": NROH,
            "called_sites": c["called_sites"],
            "called_alleles": c["called_alleles"],
            "alt_allele_count": c["alt_allele_count"],
            "hom_loci": c["hom_loci"],
            "het_loci": c["het_loci"],
            "allele_rate": safe_div(c["alt_allele_count"], c["called_alleles"]),
            "hom_loci_rate": safe_div(c["hom_loci"], c["called_sites"]),
            "het_loci_rate": safe_div(c["het_loci"], c["called_sites"]),
        }

        burden_rows.append(row)

        for comp, val in [
            ("ALT alleles", row["allele_rate"]),
            ("Homozygotes", row["hom_loci_rate"]),
            ("Heterozygotes", row["het_loci_rate"]),
        ]:
            ratio_rows.append({
                "sample": sample,
                "variant_class": cls,
                "region": region,
                "component": comp,
                "ratio": val,
                "FROH": FROH
            })

    with open(burden_file, "wt", newline="") as fh:
        fieldnames = list(burden_rows[0].keys()) if burden_rows else []
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in burden_rows:
            w.writerow(r)

    with open(ratio_file, "wt", newline="") as fh:
        fieldnames = ["sample", "variant_class", "region", "component", "ratio", "FROH"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in ratio_rows:
            w.writerow(r)

    print("[WRITE]", site_file)
    print("[WRITE]", burden_file)
    print("[WRITE]", ratio_file)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gerp4-vcf", required=True)
    ap.add_argument("--gerpneg-vcf", required=True)
    ap.add_argument("--froh", required=True)
    ap.add_argument("--roh", required=True)
    ap.add_argument("--lauto", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-samples", default="BW_01192")
    ap.add_argument("--ploidy", type=int, default=2)
    args = ap.parse_args()

    exclude = parse_exclude(args.exclude_samples)

    lgenome = read_lauto_total(args.lauto)
    froh = read_froh(args.froh)

    rohs = read_roh_segments(args.roh, exclude)
    interval_index = build_interval_index(rohs)

    counts = defaultdict(init_count)
    site_rows = []

    parse_vcf(args.gerp4_vcf, "GERP > 4", exclude, args.ploidy, interval_index, counts, site_rows)
    parse_vcf(args.gerpneg_vcf, "GERP < 0", exclude, args.ploidy, interval_index, counts, site_rows)

    write_outputs(args.out_prefix, counts, site_rows, froh, lgenome)


if __name__ == "__main__":
    main()