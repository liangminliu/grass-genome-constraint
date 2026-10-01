#!/usr/bin/env python3
"""Extract observed reference coverage from mafPairCoverage query reports."""
import argparse
import csv
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True, help="Reference MAF prefix without '*'")
    p.add_argument("--species-list", required=True, type=Path)
    p.add_argument("--results-dir", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    species = [line.split()[0] for line in args.species_list.read_text().splitlines() if line.strip()]
    if not species or len(species) != len(set(species)) or args.reference in species:
        p.error("query list must be nonempty, unique, and exclude the reference")

    rows = []
    for query in species:
        path = args.results_dir / f"{query}.tsv"
        if not path.is_file():
            p.error(f"missing coverage result: {path}")
        values = []
        for line in path.read_text().splitlines():
            fields = line.split()
            if len(fields) >= 6 and fields[0] == f"{args.reference}*":
                values.append(float(fields[5]))
        if len(values) != 1 or not 0 <= values[0] <= 1:
            p.error(f"expected one fractional ObsCoverage for {args.reference}* in {path}")
        rows.append((query, f"{100 * values[0]:.4f}"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(("Species", "ObsCoverage"))
        writer.writerows(rows)
    print(f"Wrote {len(rows)} query coverages to {args.output}")


if __name__ == "__main__":
    main()
