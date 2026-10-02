#!/usr/bin/env python3
"""Export the final bamboo and teosinte sample IDs from Supplementary Table S3."""

import argparse
from pathlib import Path

import openpyxl


PANELS = {
    "S3B_bamboo_Ho_ROH": ("selected/ped/retained_samples.txt", 63),
    "S3C_teosinte_Ho_ROH": ("selected/zmay/retained_samples.txt", 90),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("s3_workbook", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    workbook = openpyxl.load_workbook(args.s3_workbook, read_only=True, data_only=True)
    for sheet, (relative_path, expected) in PANELS.items():
        rows = workbook[sheet].iter_rows(values_only=True)
        for _ in range(3):
            next(rows)
        header = next(rows)
        iid_column = header.index("IID")
        ids = [str(row[iid_column]).strip() for row in rows if row[iid_column] is not None]
        if len(ids) != expected or len(set(ids)) != expected:
            raise ValueError(f"{sheet}: expected {expected} distinct IIDs, found {len(ids)}")
        output = args.output_dir / relative_path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(ids) + "\n", encoding="utf-8")
        print(f"{sheet}: {len(ids)} IDs -> {output}")


if __name__ == "__main__":
    main()
