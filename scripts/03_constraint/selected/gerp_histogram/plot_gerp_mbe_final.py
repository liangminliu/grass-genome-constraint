#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator, FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# -----------------------------
# Global style: MBE-like
# -----------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,

    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "xtick.color": "black",
    "ytick.color": "black",
    "text.color": "black",
    "axes.labelcolor": "black",
})


def to_superscript(n):
    """Convert integer to unicode superscript string."""
    trans = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    return str(n).translate(trans)


def style_axis(ax):
    ax.set_facecolor("white")

    # only horizontal grid
    ax.yaxis.grid(True, which="major", color="#E6E6E6", linestyle="-", linewidth=0.6)
    ax.xaxis.grid(False)

    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(0.8)

    ax.tick_params(
        axis="both",
        which="major",
        width=0.8,
        length=3.5,
        colors="black",
        pad=2
    )


def choose_unit(max_value):
    """
    Choose a scaling unit for y-axis labels.
    """
    for unit in [10**8, 10**7, 10**6, 10**5, 10**4, 10**3, 1]:
        if max_value >= unit:
            return unit
    return 1


def nice_integer_ticks(ymin, ymax, unit):
    """
    Create sparse integer ticks after scaling by unit.
    Example:
        actual values: 0 ~ 3e8
        unit: 1e8
        ticks shown: 0,1,2,3
    """
    scaled_min = ymin / unit
    scaled_max = ymax / unit
    span = scaled_max - scaled_min

    if span <= 2.5:
        step = 1
    elif span <= 5:
        step = 2
    elif span <= 9:
        step = 3
    elif span <= 15:
        step = 5
    else:
        step = max(1, int(np.ceil(span / 4.0)))

    start = np.ceil(scaled_min / step) * step
    end = np.floor(scaled_max / step) * step
    ticks_scaled = np.arange(start, end + 0.001, step)

    if ymin == 0 and (len(ticks_scaled) == 0 or ticks_scaled[0] != 0):
        ticks_scaled = np.insert(ticks_scaled, 0, 0)

    return ticks_scaled


def set_scaled_integer_yaxis(ax, ymin, ymax):
    """
    Set y-axis:
      - title: scientific notation with superscript, e.g. Frequency (×10⁸)
      - tick labels: plain numbers, e.g. 0,1,2,3
    """
    unit = choose_unit(ymax)
    ticks_scaled = nice_integer_ticks(ymin, ymax, unit)

    if len(ticks_scaled) == 0:
        ticks_scaled = np.array([ymin / unit, ymax / unit])

    ticks_actual = ticks_scaled * unit

    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_major_locator(FixedLocator(ticks_actual))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(round(y / unit))}"))

    if unit == 1:
        ylabel = "Frequency"
    else:
        exp = int(np.log10(unit))
        ylabel = f"Frequency (×10{to_superscript(exp)})"

    ax.set_ylabel(ylabel, fontsize=9)


def choose_break_point(all_hist, cds_stack, break_ratio=1.08):
    vals = np.unique(cds_stack)
    vals = vals[vals > 0]

    if len(vals) >= 3:
        ref = np.sort(vals)[-3]
    elif len(vals) == 2:
        ref = np.sort(vals)[-2]
    elif len(vals) == 1:
        ref = vals[0]
    else:
        ref = max(np.max(all_hist), np.max(cds_stack)) * 0.25

    return max(int(ref * break_ratio), 1)


def histogram_from_tsv(path, bins, bin_min, bin_max, chunksize=2_000_000, sep="\t"):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    hist = np.zeros(len(bins) - 1, dtype=np.int64)

    reader = pd.read_csv(
        path,
        sep=sep,
        header=None,
        usecols=[0],
        dtype=str,
        chunksize=chunksize,
        engine="c"
    )

    for chunk in reader:
        arr = pd.to_numeric(chunk.iloc[:, 0], errors="coerce").to_numpy()
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            continue

        arr = arr.astype(np.float32, copy=False)
        np.clip(arr, bin_min, bin_max, out=arr)
        hist += np.histogram(arr, bins=bins)[0]

    return hist


def add_panel_label(ax, label):
    if label:
        ax.text(
            -0.06, 1.02, label,
            transform=ax.transAxes,
            ha="left", va="bottom",
            fontsize=9, fontweight="bold",
            color="black"
        )


def main(args):
    fig_w = args.fig_width
    fig_h = args.fig_height

    bin_min, bin_max = -6, 6
    bins = np.arange(bin_min, bin_max + args.bin_size, args.bin_size)
    if bins[-1] < bin_max:
        bins = np.append(bins, bin_max)

    all_hist = histogram_from_tsv(
        args.all_file, bins, bin_min, bin_max,
        chunksize=args.chunksize, sep=args.sep
    )
    c_hist = histogram_from_tsv(
        args.c_cds_file, bins, bin_min, bin_max,
        chunksize=args.chunksize, sep=args.sep
    )
    d_hist = histogram_from_tsv(
        args.d_cds_file, bins, bin_min, bin_max,
        chunksize=args.chunksize, sep=args.sep
    )

    cds_stack = c_hist + d_hist
    break_point = choose_break_point(all_hist, cds_stack, break_ratio=args.break_ratio)
    upper_max = max(np.max(all_hist), np.max(cds_stack))
    upper_ylim = upper_max * 1.06

    # colors
    all_edge = "#666666"   # gray hollow bars
    c_fill =  "#464646"     #  "#9A9A9A"
    d_fill =   "#959595"      # "#D6D6D6"
    
    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1,
        sharex=True,
        figsize=(fig_w, fig_h),
        gridspec_kw={"height_ratios": [1.0, 1.0]}
    )

    for ax in (ax_top, ax_bottom):
        style_axis(ax)

    bar_width = bins[1] - bins[0]

    # -----------------------------
    # bottom panel
    # -----------------------------
    ax_bottom.bar(
        bins[:-1], c_hist,
        width=bar_width,
        align="edge",
        color=c_fill,
        edgecolor="none",
        linewidth=0.0,
        zorder=2
    )
    ax_bottom.bar(
        bins[:-1], d_hist,
        width=bar_width,
        bottom=c_hist,
        align="edge",
        color=d_fill,
        edgecolor="none",
        linewidth=0.0,
        zorder=2
    )
    ax_bottom.bar(
        bins[:-1], all_hist,
        width=bar_width,
        align="edge",
        facecolor="none",
        edgecolor=all_edge,
        linewidth=1.0,
        zorder=4
    )

    # -----------------------------
    # top panel
    # -----------------------------
    ax_top.bar(
        bins[:-1], c_hist,
        width=bar_width,
        align="edge",
        color=c_fill,
        edgecolor="none",
        linewidth=0.0,
        zorder=2
    )
    ax_top.bar(
        bins[:-1], d_hist,
        width=bar_width,
        bottom=c_hist,
        align="edge",
        color=d_fill,
        edgecolor="none",
        linewidth=0.0,
        zorder=2
    )
    ax_top.bar(
        bins[:-1], all_hist,
        width=bar_width,
        align="edge",
        facecolor="none",
        edgecolor=all_edge,
        linewidth=1.0,
        zorder=4
    )

    # y-axis: scientific notation in title, plain numbers in ticks
    set_scaled_integer_yaxis(ax_bottom, 0, break_point)
    set_scaled_integer_yaxis(ax_top, break_point, upper_ylim)

    # x-axis
    ax_bottom.set_xlim(bin_min, bin_max)
    ax_bottom.set_xticks(np.arange(-6, 7, 2))
    ax_bottom.set_xlabel("GERP score", fontsize=9)

    # hide top x labels
    ax_top.tick_params(labelbottom=False)

    # broken axis
    ax_top.spines["bottom"].set_visible(False)
    ax_bottom.spines["top"].set_visible(False)

    d = 0.012
    kwargs = dict(transform=ax_top.transAxes, color="black", clip_on=False, linewidth=0.8)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    kwargs.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    # optional break note
    if args.show_break_note:
        ax_top.text(
            0.01, 0.04, "broken y-axis",
            transform=ax_top.transAxes,
            ha="left", va="bottom",
            fontsize=7.5, color="#666666"
        )

    # overflow labels
    if args.annotate_overflow:
        ax_bottom.text(
            -6, -0.12 * break_point, "≤ -6",
            ha="center", va="top",
            fontsize=7.5, color="black", clip_on=False
        )
        ax_bottom.text(
            6, -0.12 * break_point, "≥ 6",
            ha="center", va="top",
            fontsize=7.5, color="black", clip_on=False
        )

    # legend
    handles = [
        Line2D([0], [0], color=all_edge, lw=1.0, label="All"),
        Patch(facecolor=c_fill, edgecolor="none", label="C CDS"),
        Patch(facecolor=d_fill, edgecolor="none", label="D CDS"),
    ]
    leg = ax_top.legend(
        handles=handles,
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#D0D0D0",
        framealpha=1.0,
        borderpad=0.4,
        handlelength=1.4,
        handletextpad=0.5,
        borderaxespad=0.5
    )
    leg.get_frame().set_linewidth(0.8)

    add_panel_label(ax_top, args.panel_label)

    plt.tight_layout(pad=0.5)
    plt.subplots_adjust(hspace=0.06)

    out = args.output_prefix
    fig.savefig(f"{out}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MBE-style GERP histogram with broken y-axis."
    )
    parser.add_argument("--all_file", required=True, help="All-site GERP file")
    parser.add_argument("--c_cds_file", required=True, help="C subgenome CDS GERP file")
    parser.add_argument("--d_cds_file", required=True, help="D subgenome CDS GERP file")
    parser.add_argument("--bin_size", type=float, default=0.5, help="Histogram bin width")
    parser.add_argument("--break_ratio", type=float, default=1.08, help="Multiplier for y-axis break")
    parser.add_argument("--output_prefix", type=str, default="ped_gerp_mbe_final", help="Output prefix")
    parser.add_argument("--annotate_overflow", action="store_true", help="Annotate edge bins as ≤ -6 and ≥ 6")
    parser.add_argument("--show_break_note", action="store_true", help="Add a small 'broken y-axis' note")
    parser.add_argument("--panel_label", type=str, default="", help="Panel label, e.g. A or B")
    parser.add_argument("--chunksize", type=int, default=2000000, help="Rows per chunk")
    parser.add_argument("--sep", type=str, default="\t", help="Field separator, default is tab")
    parser.add_argument("--fig_width", type=float, default=4.8, help="Figure width in inches")
    parser.add_argument("--fig_height", type=float, default=4.2, help="Figure height in inches")
    args = parser.parse_args()
    main(args)