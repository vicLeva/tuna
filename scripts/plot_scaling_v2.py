#!/usr/bin/env python3
"""n-files scaling — tuna/KMC/FastK (old) + tuna_rob (rob_v2).

Reads:
  benchmark/results/bench_scaling_*/bench_scaling.csv  — most recent 3-tool CSV
  benchmark/results/rob_v2/bench_scaling.csv           — tuna_rob scaling CSV

Output PNGs (dpi=500) in benchmark/results/rob_v2/:
  scaling_ecoli.png  scaling_human.png  scaling_combined.png
"""

import csv, sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import ConnectionPatch

RESULTS = Path(__file__).parent.parent / "benchmark" / "results"
ROB_CSV = RESULTS / "rob_v2" / "bench_scaling.csv"
OUT_DIR = RESULTS / "rob_v2"

TOOLS = ["tuna", "kmc", "fastk", "tuna_rob"]
STYLE = {
    "tuna":     dict(color="#2166ac", marker="o",  label="tuna (dev)"),
    "kmc":      dict(color="#d6604d", marker="s",  label="KMC3"),
    "fastk":    dict(color="#4dac26", marker="^",  label="FastK"),
    "tuna_rob": dict(color="#e6ab02", marker="v",  label="tuna_rob", linewidth=2.2),
}

DATASET_META = {
    "ecoli": dict(title="E. coli  (k=31, assemblies)"),
    "human": dict(title="Human    (k=31, assemblies)"),
}


# ── formatters ────────────────────────────────────────────────────────────────

def _fmt_time(v, _):
    if v <= 0:  return "0 s"
    if v < 10:  return f"{v:.1f} s"
    return f"{v:.0f} s"

def _fmt_n(v, _):
    v = int(v)
    return f"{v//1000}k" if v >= 1000 else str(v)


# ── loading ───────────────────────────────────────────────────────────────────

def load_old(csv_path):
    """data[ds][tool] = sorted [(n_files, wall_s, rss_mb)]."""
    data = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, tool, n = row["dataset"], row["tool"], int(row["n_files"])
            wall = float(row["wall_s"])   if row.get("wall_s","")  not in ("","na") else None
            rss  = float(row["rss_mb"])   if row.get("rss_mb","")  not in ("","na") else None
            if wall is not None:
                data[ds][tool].append((n, wall, rss))
    for ds in data:
        for t in data[ds]:
            data[ds][t].sort()
    return data


def load_rob(csv_path):
    """Returns data in the same shape, tool='tuna_rob'."""
    data = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, n = row["dataset"], int(row["n_files"])
            wall = float(row["wall_s"])  if row.get("wall_s","")  not in ("","na") else None
            rss  = float(row["rss_mb"])  if row.get("rss_mb","")  not in ("","na") else None
            if wall is not None:
                data[ds]["tuna_rob"].append((n, wall, rss))
    for ds in data:
        data[ds]["tuna_rob"].sort()
    return data


def merge(old, new):
    merged = defaultdict(lambda: defaultdict(list))
    for src in (old, new):
        for ds in src:
            for tool in src[ds]:
                merged[ds][tool] = src[ds][tool]
    return merged


def pts(data, ds, tool, field):
    result = []
    for n, wall, rss in data[ds].get(tool, []):
        v = wall if field == "wall" else rss
        if v is not None:
            result.append((n, v))
    return zip(*result) if result else ([], [])


# ── axis helpers ──────────────────────────────────────────────────────────────

def _style_xaxis(ax, xs_all):
    ax.set_xscale("log")
    ticks = sorted(set(int(x) for x in xs_all))
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_fmt_n))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.tick_params(axis="x", labelsize=7, rotation=45)


def draw_time_ax(ax, data, dataset, add_legend=False):
    all_xs = []
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = pts(data, dataset, tool, "wall")
        xs, ys = list(xs), list(ys)
        if not xs:
            continue
        all_xs.extend(xs)
        lw = st.get("linewidth", 1.8)
        ax.plot(xs, ys, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=lw, markersize=6, zorder=3)
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
        xs, ys = pts(data, "ecoli", tool, "wall")
        pairs = [(x, y) for x, y in zip(xs, ys) if x <= 50]
        if not pairs:
            continue
        xi, yi = zip(*pairs)
        lw = st.get("linewidth", 1.4)
        axins.plot(xi, yi, color=st["color"], marker=st["marker"],
                   linewidth=lw, markersize=4, zorder=3)
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
            xyB=(xfrac, 0.0),   coordsB="axes fraction", axesB=axins,
            color="grey", alpha=0.55, linewidth=0.8, zorder=5)
        ax.get_figure().add_artist(conn)


def draw_rss_ax(ax, data, dataset, add_legend=False):
    all_xs = []
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = pts(data, dataset, tool, "rss")
        xs, ys = list(xs), list(ys)
        if not xs:
            continue
        all_xs.extend(xs)
        ys_gb = [v / 1024 for v in ys]
        lw = st.get("linewidth", 1.8)
        ax.plot(xs, ys_gb, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=lw, markersize=6, zorder=3)
    if all_xs:
        _style_xaxis(ax, all_xs)
    ax.set_xlabel("Number of input files", fontsize=9)
    ax.set_ylabel("Peak RSS (GB)", fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper left")


# ── figures ───────────────────────────────────────────────────────────────────

def plot_dataset(data, dataset, outprefix):
    meta = DATASET_META[dataset]
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    fig.suptitle(meta["title"], fontsize=11, fontweight="bold")
    draw_time_ax(ax_t, data, dataset, add_legend=True)
    if dataset == "ecoli":
        add_ecoli_zoom(ax_t, data)
    draw_rss_ax(ax_r, data, dataset)
    out = f"{outprefix}_{dataset}.png"
    fig.savefig(out, dpi=500)
    plt.close(fig)
    print(f"  saved  {out}")


def plot_combined(data, outprefix):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting — scaling with number of input files",
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

    out = f"{outprefix}_combined.png"
    fig.savefig(out, dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    candidates = sorted(RESULTS.glob("bench_scaling_*/bench_scaling.csv"))
    if not candidates:
        sys.exit("error: no bench_scaling_*/bench_scaling.csv under benchmark/results/")
    old_csv = candidates[-1]

    print(f"Old CSV : {old_csv}")
    print(f"Rob CSV : {ROB_CSV}")
    print(f"Output  : {OUT_DIR}/scaling_*.png")

    data = merge(load_old(old_csv), load_rob(ROB_CSV))

    outprefix = str(OUT_DIR / "scaling")
    plot_combined(data, outprefix)


if __name__ == "__main__":
    main()
