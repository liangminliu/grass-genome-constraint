#!/usr/bin/env python3
"""Calculate cumulative alignment depth from one reference's chromosome BEDs.

BED column 6 must contain the number of aligned species (k). Each BED file has
one header row. Default weights are interval lengths in bp, matching the
manuscript's site-level denominator. Use --weight rows to reproduce the
historical all.py record-count calculation.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed-dir", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--weight", choices=("bases", "rows"), default="bases")
    p.add_argument("--max-k", type=int, default=None,
                   help="Optional common maximum k for three-species plots")
    args = p.parse_args()
    files = sorted(args.bed_dir.glob("*.bed"))
    if not files:
        p.error(f"no .bed files in {args.bed_dir}")

    counts = defaultdict(int)
    n_rows = n_bases = n_nonunit = 0
    for path in files:
        with path.open(encoding="utf-8") as handle:
            next(handle, None)  # historical input format has one header row
            for line_number, line in enumerate(handle, start=2):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    p.error(f"{path}:{line_number}: expected at least 6 BED columns")
                try:
                    start, end, k = int(fields[1]), int(fields[2]), int(fields[5])
                except ValueError:
                    p.error(f"{path}:{line_number}: invalid start, end, or k")
                if start < 0 or end <= start or k < 0:
                    p.error(f"{path}:{line_number}: invalid interval or k")
                length = end - start
                counts[k] += length if args.weight == "bases" else 1
                n_rows += 1
                n_bases += length
                n_nonunit += length != 1
    if not n_rows or max(counts) < 1:
        p.error("no positive alignment-depth values found")
    if args.max_k is not None and args.max_k < max(counts):
        p.error(f"--max-k {args.max_k} is below observed k={max(counts)}")

    max_k = args.max_k or max(counts)
    denominator = sum(counts.values())
    surviving = denominator
    values = []
    for k in range(1, max_k + 1):
        surviving -= counts.get(k - 1, 0)
        values.append((k, f"{100 * surviving / denominator:.8f}"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(("k", "All"))
        writer.writerows(values)
    print(f"Wrote {args.output}; BED rows={n_rows}, interval bp={n_bases}, "
          f"non-unit intervals={n_nonunit}, denominator={denominator} {args.weight}")


if __name__ == "__main__":
    main()
