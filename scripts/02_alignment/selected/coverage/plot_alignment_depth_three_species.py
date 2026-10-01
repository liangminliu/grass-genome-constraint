#!/usr/bin/env python3
"""Plot three cumulative alignment-depth curves for manuscript Fig. 2B."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_curve(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any("k" not in row or "All" not in row for row in rows):
        raise ValueError(f"expected k/All columns in {path}")
    x = [int(row["k"]) for row in rows]
    y = [float(row["All"]) for row in rows]
    if x != list(range(1, len(x) + 1)) or any(value < 0 or value > 100 for value in y):
        raise ValueError(f"invalid cumulative curve in {path}")
    if any(y[n] > y[n-1] + 1e-6 for n in range(1, len(y))):
        raise ValueError(f"curve increases with k in {path}")
    return x, y


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ped", required=True, type=Path)
    p.add_argument("--zmay", required=True, type=Path)
    p.add_argument("--atau", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    inputs = (
        ("P. edulis", args.ped, "#404040", "o"),
        ("Z. mays", args.zmay, "#6A6A6A", "s"),
        ("A. tauschii", args.atau, "#9A9A9A", "^"),
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    maximum = 0
    for name, path, color, marker in inputs:
        x, y = read_curve(path)
        maximum = max(maximum, x[-1])
        ax.plot(x, y, linewidth=1, marker=marker, markersize=3,
                color=color, label=name, alpha=0.9)
    ax.set_xlabel("Minimum number of aligned species (k)")
    ax.set_ylabel("Cumulative alignment proportion (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(1, maximum)
    ax.set_xticks(list(range(10, maximum + 1, 10)))
    ax.grid(True, linestyle=":", alpha=0.3)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
