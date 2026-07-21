#!/usr/bin/env python3
"""Plot n-files scaling sweep from bench_scaling.sh output.

Usage:
    python plot_scaling.py <results_dir> [out.png]

<results_dir> is the directory produced by bench_scaling.sh (e.g. after
scp'ing it back from the cluster) — must contain bench_scaling.csv.
"""

import sys
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import ConnectionPatch

# ── per-tool visual identity (matches plot_all_tools.py) ──────────────────

STYLE = {
    "tuna":  dict(color="#2166ac", marker="o", label="tuna"),
    "kmc":   dict(color="#d6604d", marker="s", label="KMC3"),
    "fastk": dict(color="#4dac26", marker="^", label="FastK"),
}
TOOLS = ["tuna", "kmc", "fastk"]

DATASET_META = {
    "ecoli": dict(title="E. coli  (k=31, m=21, assemblies)"),
    "human": dict(title="Human    (k=31, m=21, assemblies)"),
}

# ── tick formatters ───────────────────────────────────────────────────────────


def _fmt_time(v, _):
    if v <= 0:
        return "0 s"
    if v < 10:
        return f"{v:.1f} s"
    return f"{v:.0f} s"


def _fmt_n(v, _):
    v = int(v)
    return f"{v//1000}k" if v >= 1000 else str(v)


# ── data loading ──────────────────────────────────────────────────────────────


def load(csv_path):
    """data[dataset][tool] = sorted list of (n_files, wall_s, rss_mb)."""
    data = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, tool, n = row["dataset"], row["tool"], int(row["n_files"])

            def _f(k):
                v = row.get(k, "")
                return float(v) if v not in ("", "na", "timeout") else None

            wall = _f("wall_s")
            rss = _f("rss_mb")
            if wall is not None:
                data[ds][tool].append((n, wall, rss))
    for ds in data:
        for tool in data[ds]:
            data[ds][tool].sort()
    return data


def xs_ys(data, ds, tool, field):
    fi = 0 if field == "wall" else 1
    pts = [(n, wall, rss) for n, wall, rss in data[ds].get(tool, [])
           if (wall if field == "wall" else rss) is not None]
    if not pts:
        return [], []
    xs = [p[0] for p in pts]
    ys = [p[fi + 1] for p in pts]
    return xs, ys


# ── axis drawing ──────────────────────────────────────────────────────────────


def _style_xaxis(ax, xs_all):
    ax.set_xscale("log")
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
        if not xs:
            continue
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
    axins = ax.inset_axes([0.08, 0.45, 0.45, 0.52])
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(data, "ecoli", tool, "wall")
        pairs = [(x, y) for x, y in zip(xs, ys) if x <= 50]
        if not pairs:
            continue
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

    indicator = ax.indicate_inset_zoom(axins, edgecolor="grey", alpha=0.55)
    for c in indicator.connectors:
        c.set_visible(False)

    x0, x1 = axins.get_xlim()
    y_top = axins.get_ylim()[1]
    for xdata, xfrac in [(x0, 0.0), (x1, 1.0)]:
        conn = ConnectionPatch(
            xyA=(xdata, y_top), coordsA="data", axesA=ax,
            xyB=(xfrac, 0.0), coordsB="axes fraction", axesB=axins,
            color="grey", alpha=0.55, linewidth=0.8, zorder=5,
        )
        ax.get_figure().add_artist(conn)


def draw_rss_ax(ax, data, dataset, add_legend=False):
    all_xs = []
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(data, dataset, tool, "rss")
        if not xs:
            continue
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


# ── combined figure ───────────────────────────────────────────────────────────


def plot_combined(data, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting — scaling with number of input files  (k=31, m=21, count-only)",
                 fontsize=11, fontweight="bold")

    row_labels = {"ecoli": "E. coli", "human": "Human"}
    for row, ds in enumerate(["ecoli", "human"]):
        draw_time_ax(axes[row, 0], data, ds, add_legend=False)
        if ds == "ecoli":
            add_ecoli_zoom(axes[row, 0], data)
        draw_rss_ax(axes[row, 1], data, ds, add_legend=False)
        axes[row, 0].set_ylabel(f"{row_labels[ds]}\nWall-clock time", fontsize=9)
        axes[row, 1].set_ylabel("Peak RSS (GB)", fontsize=9)

    seen, handles, labels = set(), [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in seen:
                handles.append(h); labels.append(l); seen.add(l)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.13, top=0.92)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.01),
               ncol=len(handles), fontsize=8.5, framealpha=0.85)

    fig.savefig(out_path, dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: plot_scaling.py <results_dir> [out.png]")

    results_dir = Path(sys.argv[1])
    csv_path = results_dir / "bench_scaling.csv"
    if not csv_path.is_file():
        sys.exit(f"error: not found: {csv_path}")

    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else results_dir / "scaling_combined.png"

    print(f"CSV:    {csv_path}")
    print(f"Output: {out_path}")

    data = load(csv_path)
    plot_combined(data, out_path)


if __name__ == "__main__":
    main()
