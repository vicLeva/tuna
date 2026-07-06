#!/usr/bin/env python3
"""Plot n-files scaling sweep from bench_scaling.sh output.

Produces figures matching the style of plot_all_tools.py:
  <prefix>_ecoli.png   — wall time + RSS for E. coli
  <prefix>_human.png   — wall time + RSS for Human
  <prefix>_combined.png — 2×2 combined figure

Usage:
    python scripts/plot_scaling.py [bench_scaling.csv] [--out PREFIX]
"""

import sys
import csv
import argparse
import math
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── per-tool visual identity (matches plot_all_tools.py exactly) ──────────────

STYLE = {
    "tuna":  dict(color="#2166ac", marker="o",  label="tuna"),
    "kmc":   dict(color="#d6604d", marker="s",  label="KMC3"),
    "fastk": dict(color="#4dac26", marker="^",  label="FastK"),
}
TOOLS = ["tuna", "kmc", "fastk"]

DATASET_META = {
    "ecoli": dict(title="E. coli  (k=31, assemblies)"),
    "human": dict(title="Human    (k=31, assemblies)"),
}

# ── tick formatters (matches plot_all_tools.py) ───────────────────────────────

def _fmt_time(v, _):
    if v <= 0:  return "0 s"
    if v < 10:  return f"{v:.1f} s"
    return f"{v:.0f} s"

def _fmt_n(v, _):
    """Compact label for n_files log axis."""
    v = int(v)
    if v >= 1000: return f"{v//1000}k"
    return str(v)

# ── data loading ──────────────────────────────────────────────────────────────

def load(csv_path):
    """data[dataset][tool] = sorted list of (n_files, wall_s, rss_mb)."""
    data = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, tool, n = row["dataset"], row["tool"], int(row["n_files"])
            def _f(k):
                v = row.get(k, "")
                return float(v) if v not in ("", "na") else None
            wall = _f("wall_s")
            rss  = _f("rss_mb")
            if wall is not None:
                data[ds][tool].append((n, wall, rss))
    for ds in data:
        for tool in data[ds]:
            data[ds][tool].sort()
    return data

def xs_ys(data, ds, tool, field):
    """field: 'wall' or 'rss'. Returns (xs, ys) with None values dropped."""
    fi = 0 if field == "wall" else 1
    pts = [(n, wall, rss) for n, wall, rss in data[ds].get(tool, [])
           if (wall if field == "wall" else rss) is not None]
    if not pts: return [], []
    xs = [p[0] for p in pts]
    ys = [p[fi + 1] for p in pts]
    return xs, ys

# ── axis drawing ──────────────────────────────────────────────────────────────

def _style_xaxis(ax, xs_all):
    """Log x-axis with clean integer tick labels."""
    ax.set_xscale("log")
    # pick explicit ticks from the data itself
    ticks = sorted(set(xs_all))
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_fmt_n))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.tick_params(axis="x", labelsize=7, rotation=45)

def draw_time_ax(ax, data, dataset, add_legend=False):
    all_xs = []
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(data, dataset, tool, "wall")
        if not xs: continue
        all_xs.extend(xs)
        ax.plot(xs, ys, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=1.8, markersize=6, zorder=3)
    if all_xs:
        _style_xaxis(ax, all_xs)
    ax.set_xlabel("Number of input files", fontsize=9)
    ax.set_ylabel("Wall-clock time", fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper left")

def add_ecoli_zoom(ax, data):
    """Inset zooming into n=1..50 on the E. coli wall-time panel."""
    from matplotlib.patches import ConnectionPatch

    axins = ax.inset_axes([0.08, 0.45, 0.45, 0.52])
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(data, "ecoli", tool, "wall")
        pairs = [(x, y) for x, y in zip(xs, ys) if x <= 50]
        if not pairs: continue
        xi, yi = zip(*pairs)
        axins.plot(xi, yi, color=st["color"], marker=st["marker"],
                   linewidth=1.4, markersize=4, zorder=3)
    axins.set_xscale("log")
    axins.set_xticks([1, 2, 5, 10, 20, 50])
    axins.xaxis.set_major_formatter(mticker.FuncFormatter(_fmt_n))
    axins.xaxis.set_minor_locator(mticker.NullLocator())
    axins.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    axins.tick_params(labelsize=6)
    axins.grid(axis="y", linestyle=":", alpha=0.35)
    axins.spines[["top", "right"]].set_visible(False)

    # draw the indicator rectangle but hide the default connectors
    indicator = ax.indicate_inset_zoom(axins, edgecolor="grey", alpha=0.55)
    for c in indicator.connectors:
        c.set_visible(False)

    # custom connectors: top-left and top-right of indicator → bottom corners of inset
    x0, x1 = axins.get_xlim()
    y_top   = axins.get_ylim()[1]   # top of the zoomed data range = top of indicator box
    for xdata, xfrac in [(x0, 0.0), (x1, 1.0)]:
        conn = ConnectionPatch(
            xyA=(xdata, y_top), coordsA="data", axesA=ax,
            xyB=(xfrac, 0.0),   coordsB="axes fraction", axesB=axins,
            color="grey", alpha=0.55, linewidth=0.8, zorder=5,
        )
        ax.get_figure().add_artist(conn)

def draw_rss_ax(ax, data, dataset, add_legend=False):
    all_xs = []
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(data, dataset, tool, "rss")
        if not xs: continue
        all_xs.extend(xs)
        ys_gb = [v / 1024 for v in ys]
        ax.plot(xs, ys_gb, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=1.8, markersize=6, zorder=3)
    if all_xs:
        _style_xaxis(ax, all_xs)
    ax.set_xlabel("Number of input files", fontsize=9)
    ax.set_ylabel("Peak RSS (GB)", fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper left")

# ── per-dataset figure ────────────────────────────────────────────────────────

def plot_dataset(data, dataset, outprefix):
    meta = DATASET_META[dataset]
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9, 3.8),
                                      constrained_layout=True)
    fig.suptitle(meta["title"], fontsize=11, fontweight="bold")
    draw_time_ax(ax_t, data, dataset, add_legend=True)
    if dataset == "ecoli":
        add_ecoli_zoom(ax_t, data)
    draw_rss_ax(ax_r, data, dataset, add_legend=False)
    out = f"{outprefix}_{dataset}.png"
    fig.savefig(out, dpi=500)
    plt.close(fig)
    print(f"  saved  {out}")

# ── combined 2×2 figure ───────────────────────────────────────────────────────

def plot_combined(data, outprefix):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting — scaling with number of input files  (k=31, count + text dump)",
                 fontsize=11, fontweight="bold")

    row_labels = {"ecoli": "E. coli", "human": "Human"}

    for row, ds in enumerate(["ecoli", "human"]):
        draw_time_ax(axes[row, 0], data, ds, add_legend=False)
        if ds == "ecoli":
            add_ecoli_zoom(axes[row, 0], data)
        draw_rss_ax(axes[row, 1], data, ds, add_legend=False)
        axes[row, 0].set_ylabel(f"{row_labels[ds]}\nWall-clock time", fontsize=9)
        axes[row, 1].set_ylabel("Peak RSS (GB)", fontsize=9)

    # shared legend below
    seen, handles, labels = set(), [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in seen:
                handles.append(h); labels.append(l); seen.add(l)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.13, top=0.92)
    fig.legend(handles, labels,
               loc="lower center", bbox_to_anchor=(0.5, 0.01),
               ncol=len(handles), fontsize=8.5, framealpha=0.85)

    out = f"{outprefix}_combined.png"
    fig.savefig(out, dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  {out}")

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", nargs="?", type=Path,
                    help="bench_scaling.csv (default: most recent in benchmark/results/)")
    ap.add_argument("--out", default=None,
                    help="output prefix (default: same dir as csv, stem=scaling)")
    args = ap.parse_args()

    csv_path = args.csv
    if csv_path is None:
        candidates = sorted(
            Path(__file__).parent.parent.glob(
                "benchmark/results/bench_scaling_*/bench_scaling.csv"))
        if not candidates:
            sys.exit("error: no bench_scaling.csv found under benchmark/results/")
        csv_path = candidates[-1]

    outprefix = args.out or str(csv_path.parent / "scaling")

    print(f"CSV:    {csv_path}")
    print(f"Output: {outprefix}_*.png")

    data = load(csv_path)

    for ds in ["ecoli", "human"]:
        if ds in data:
            plot_dataset(data, ds, outprefix)

    plot_combined(data, outprefix)

if __name__ == "__main__":
    main()
