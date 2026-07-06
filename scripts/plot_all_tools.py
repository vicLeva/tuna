#!/usr/bin/env python3
"""
Parse /usr/bin/time -v benchmark results from bench_default.sh and produce
publication-quality plots (wall time + peak RSS vs. thread count).

Usage:
    python plot_bench_default.py results.logs [--out PREFIX]

Outputs  (default prefix = "bench"):
    bench_ecoli.png
    bench_human.png
    bench_combined.png   (2×2 grid, all datasets × metrics)
"""

import re
import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


# ── per-tool visual identity ────────────────────────────────────────────────
TOOLS = ["tuna", "kmc", "fastk", "kfc"]

STYLE = {
    "tuna":  dict(color="#2166ac", marker="o",  label="tuna"),
    "kmc":   dict(color="#d6604d", marker="s",  label="KMC3"),
    "fastk": dict(color="#4dac26", marker="^",  label="FastK"),
    "kfc":   dict(color="#8856a7", marker="D",  label="KFC"),
}

THREADS = [1, 2, 4, 6, 8, 10, 16, 22, 28, 32]


# ── parsing ─────────────────────────────────────────────────────────────────

def _wall_to_sec(s: str) -> float:
    """'h:mm:ss.cc'  |  'm:ss.cc'  →  float seconds."""
    parts = s.strip().split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def parse_results(path: Path) -> list[dict]:
    """Return one dict per run block found in *path*."""
    hdr_re  = re.compile(
        r"=== TOOL=(\w+)\s+DS=(\w+)\s+FILE=\S+\s+T=(\d+)\s+DATE=")
    wall_re = re.compile(r"Elapsed \(wall clock\) time.*?:\s+([\d:\.]+)")
    rss_re  = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
    exit_re = re.compile(r"Exit status:\s+(\d+)")

    text   = path.read_text(errors="replace")
    chunks = re.split(r"(?=^=== TOOL=)", text, flags=re.MULTILINE)

    records = []
    for chunk in chunks:
        m = hdr_re.search(chunk)
        if not m:
            continue
        tool, ds, threads = m.group(1).lower(), m.group(2).lower(), int(m.group(3))

        wm = wall_re.search(chunk)
        rm = rss_re.search(chunk)
        em = exit_re.search(chunk)

        records.append(dict(
            tool    = tool,
            dataset = ds,
            threads = threads,
            wall_s  = _wall_to_sec(wm.group(1)) if wm else None,
            rss_kb  = int(rm.group(1))           if rm else None,
            ok      = em is not None and int(em.group(1)) == 0,
        ))
    return records


def lookup(records, tool, dataset, threads, field):
    """Return (value, ok) or (None, None) when not found."""
    hits = [r for r in records
            if r["tool"] == tool and r["dataset"] == dataset
            and r["threads"] == threads]
    if not hits:
        return None, None
    r = hits[0]
    return r.get(field), r["ok"]


# ── plot helpers ─────────────────────────────────────────────────────────────

def _fmt_time(v, _):
    """Tick formatter: always in seconds."""
    if v <= 0:
        return "0 s"
    if v < 10:
        return f"{v:.1f} s"
    return f"{v:.0f} s"


def _fmt_time_sparse(v, pos):
    """Log-axis formatter: gridlines at every sub-decade, labels only at 1×/2×/5×."""
    if v <= 0:
        return ""
    import math
    exp = math.floor(math.log10(v))
    m = v / 10 ** exp          # mantissa in [1, 10)
    if any(abs(m - x) < 0.05 for x in (1, 2, 5)):
        return _fmt_time(v, pos)
    return ""


def _fmt_rss(v, _, unit_gb=True):
    if unit_gb:
        return f"{v/1024**2:.0f}" if v >= 1024**2 else f"{v/1024:.0f} MB"
    return f"{v/1024:.0f}"


def draw_time_ax(ax, records, dataset, logscale=False, add_legend=True):
    """Line plot of wall time vs threads for one dataset."""
    any_crash = False

    for tool in TOOLS:
        st = STYLE[tool]
        xs_ok, ys_ok   = [], []
        xs_bad, ys_bad = [], []

        for t in THREADS:
            val, ok = lookup(records, tool, dataset, t, "wall_s")
            if val is None:
                continue
            if ok:
                xs_ok.append(t);  ys_ok.append(val)
            else:
                xs_bad.append(t); ys_bad.append(val)
                any_crash = True

        if xs_ok:
            ax.plot(xs_ok, ys_ok,
                    color=st["color"], marker=st["marker"],
                    label=st["label"],
                    linewidth=1.8, markersize=6, zorder=3)
        if xs_bad:
            ax.plot(xs_bad, ys_bad,
                    color=st["color"], marker="x",
                    linestyle="--", linewidth=1.2,
                    markersize=9, markeredgewidth=2.0, zorder=3,
                    label=f"{st['label']} (crashed)")
            # connect the dashed segment to the last valid ok point
            if xs_ok and xs_bad[0] > xs_ok[-1]:
                ax.plot([xs_ok[-1], xs_bad[0]],
                        [ys_ok[-1], ys_bad[0]],
                        color=st["color"], linestyle="--",
                        linewidth=1.2, zorder=2)

    ax.set_xlabel("Threads", fontsize=9)
    ylabel = "Wall-clock time (log scale)" if logscale else "Wall-clock time"
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xticks(THREADS)
    ax.set_xscale("linear")
    ax.xaxis.set_minor_locator(mticker.NullLocator())

    if logscale:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(
            mticker.LogLocator(base=10, subs=[1, 2, 3, 4, 5, 6, 7, 8, 9], numticks=20))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time_sparse))
        ax.yaxis.set_minor_locator(mticker.NullLocator())
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_ylim(bottom=90)   # ensure the 100 s gridline/label is visible
    else:
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))

    if any_crash:
        ax.annotate("✕ = crash (Tabex SIGSEGV)", xy=(0.97, 0.03),
                    xycoords="axes fraction", ha="right", va="bottom",
                    fontsize=7, color="grey", style="italic")

    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)

    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper right")


def draw_rss_ax(ax, records, dataset, unit_gb=True, add_legend=False):
    """Line plot of peak RSS vs threads for one dataset."""
    unit_div  = 1024**2 if unit_gb else 1024      # kB → GB  or  kB → MB
    unit_label = "GB" if unit_gb else "MB"

    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = [], []

        for t in THREADS:
            val, ok = lookup(records, tool, dataset, t, "rss_kb")
            if val is None:
                continue
            xs.append(t)
            ys.append(val / unit_div)

        if xs:
            ax.plot(xs, ys,
                    color=st["color"], marker=st["marker"],
                    label=st["label"],
                    linewidth=1.8, markersize=6, zorder=3)

    ax.set_xlabel("Threads", fontsize=9)
    ax.set_ylabel(f"Peak RSS ({unit_label})", fontsize=9)
    ax.set_xticks(THREADS)
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)

    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper right")


# ── per-dataset figure ───────────────────────────────────────────────────────

DATASET_META = {
    "ecoli": dict(title="E. coli (k=31, 1 genome)",
                  logtime=False, rss_gb=True),
    "human": dict(title="Human (k=31, 1 genome)",
                  logtime=True,  rss_gb=True),
}


def plot_dataset(records, dataset, meta, outprefix):
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9, 3.8),
                                      constrained_layout=True)
    fig.suptitle(meta["title"], fontsize=11, fontweight="bold")

    draw_time_ax(ax_t, records, dataset,
                 logscale=meta["logtime"], add_legend=True)
    draw_rss_ax(ax_r, records, dataset,
                unit_gb=meta["rss_gb"], add_legend=False)

    fig.savefig(f"{outprefix}_{dataset}.png", dpi=500)
    plt.close(fig)
    print(f"  saved  {outprefix}_{dataset}.png ")


# ── combined 2×2 figure ──────────────────────────────────────────────────────

def plot_combined(records, outprefix):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting benchmark  (k = 31, count + text dump)",
                 fontsize=11, fontweight="bold")

    datasets   = ["ecoli", "human"]
    row_labels = ["E. coli", "Human"]

    for row, (ds, rl) in enumerate(zip(datasets, row_labels)):
        meta = DATASET_META[ds]

        draw_time_ax(axes[row, 0], records, ds,
                     logscale=False,
                     add_legend=False)
        draw_rss_ax(axes[row, 1], records, ds,
                    unit_gb=meta["rss_gb"],
                    add_legend=False)

        axes[row, 0].set_ylabel(f"{rl}\nWall-clock time", fontsize=9)
        axes[row, 1].set_ylabel(
            f"Peak RSS ({'GB' if meta['rss_gb'] else 'MB'})", fontsize=9)

    # ── zoom inset on human wall-time: [100 s, 500 s] detail ─────────────────
    ax_ht = axes[1, 0]
    axins = ax_ht.inset_axes([0.28, 0.24, 0.69, 0.56])
    for tool in ["tuna", "kmc", "fastk"]:
        st = STYLE[tool]
        xs, ys = [], []
        for t in THREADS:
            val, ok = lookup(records, tool, "human", t, "wall_s")
            if val is not None and ok:
                xs.append(t); ys.append(val)
        if xs:
            axins.plot(xs, ys, color=st["color"], marker=st["marker"],
                       linewidth=1.4, markersize=4, zorder=3)
    axins.set_xlim(THREADS[0] - 0.5, THREADS[-1] + 0.5)
    axins.set_ylim(90, 550)
    axins.set_xticks([1, 4, 8, 16, 22, 28])
    axins.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    axins.tick_params(labelsize=6)
    axins.grid(axis="y", linestyle=":", alpha=0.35)
    axins.spines[["top", "right"]].set_visible(False)
    ax_ht.indicate_inset_zoom(axins, edgecolor="grey", alpha=0.55)

    # ── single shared legend below all subplots ───────────────────────────────
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

    fig.savefig(f"{outprefix}_combined.png", dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  {outprefix}_combined.png ")


# ── per-thread comparison figure ─────────────────────────────────────────────

def plot_per_thread(records, outprefix):
    """One figure per thread count: grouped bar chart, both datasets side-by-side."""
    datasets = ["ecoli", "human"]
    n_tools  = len(TOOLS)

    for t in THREADS:
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8),
                                  constrained_layout=True)
        fig.suptitle(f"k-mer counting benchmark  –  {t} thread{'s' if t>1 else ''}",
                     fontsize=11, fontweight="bold")

        for col, ds in enumerate(datasets):
            ax   = axes[col]
            meta = DATASET_META[ds]
            unit_div  = 1024**2 if meta["rss_gb"] else 1024
            unit_label = "GB" if meta["rss_gb"] else "MB"

            x     = np.arange(n_tools)
            width = 0.35

            times, rsss, colors, hatches, tick_labels = [], [], [], [], []

            for tool in TOOLS:
                st   = STYLE[tool]
                wall, ok  = lookup(records, tool, ds, t, "wall_s")
                rss,  ok2 = lookup(records, tool, ds, t, "rss_kb")
                times.append(wall if wall is not None else 0)
                rsss.append(rss / unit_div if rss is not None else 0)
                colors.append(st["color"])
                hatches.append("///" if (not ok and wall is not None) else "")
                tick_labels.append(st["label"])

            # time bars
            bars_t = ax.bar(x - width/2, times, width,
                            color=colors,
                            hatch=[h or None for h in hatches],
                            edgecolor="white", linewidth=0.6,
                            label="Wall time")

            ax2 = ax.twinx()
            bars_r = ax2.bar(x + width/2, rsss, width,
                             color=colors, alpha=0.4,
                             edgecolor="white", linewidth=0.6,
                             label="Peak RSS")

            ax.set_xticks(x)
            ax.set_xticklabels(tick_labels, fontsize=9)
            ax.set_ylabel("Wall-clock time", fontsize=9)
            ax2.set_ylabel(f"Peak RSS ({unit_label})", fontsize=9)

            ds_name = "E. coli" if ds == "ecoli" else "Human"
            ax.set_title(ds_name, fontsize=10)

            if meta["logtime"]:
                ax.set_yscale("log")
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
            ax.grid(axis="y", linestyle=":", alpha=0.45)
            ax.spines[["top"]].set_visible(False)

            # legend patches
            from matplotlib.patches import Patch
            legend_elems = [Patch(fc=STYLE[tool]["color"], label=STYLE[tool]["label"])
                            for tool in TOOLS]
            legend_elems += [Patch(fc="grey", hatch="///", label="crashed")]
            ax.legend(handles=legend_elems, fontsize=7,
                      framealpha=0.85, loc="upper right")

        fig.savefig(f"{outprefix}_t{t}.png", dpi=500)
        plt.close(fig)
        print(f"  saved  {outprefix}_t{t}.png ")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results", type=Path,
                    help="results.logs produced by bench_default.sh")
    ap.add_argument("--out", default="bench",
                    help="output file prefix (default: bench)")
    ap.add_argument("--per-thread", action="store_true",
                    help="also produce one grouped-bar figure per thread count")
    args = ap.parse_args()

    if not args.results.exists():
        sys.exit(f"error: {args.results} not found")

    records = parse_results(args.results)
    if not records:
        sys.exit("error: no run blocks found in the log file")

    print(f"Parsed {len(records)} run(s) from {args.results}")

    outprefix = args.out

    # per-dataset line plots (main output)
    for ds, meta in DATASET_META.items():
        plot_dataset(records, ds, meta, outprefix)

    # combined 2×2
    plot_combined(records, outprefix)

    # optional: one bar chart per thread count
    if args.per_thread:
        plot_per_thread(records, outprefix)


if __name__ == "__main__":
    main()
