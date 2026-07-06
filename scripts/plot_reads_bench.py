#!/usr/bin/env python3
"""
Plot results from bench_all_tools_reads.sh (fixed T=8, FASTQ reads).

Usage:
    python scripts/plot_reads_bench.py  <results.log|dir> [<results.log|dir> ...]
                                        [--label L1 L2 ...]
                                        [--out PREFIX]

Multiple log files are overlaid as separate bar groups (e.g. local vs cluster).
If a single log is given the label defaults to the parent directory name.

Outputs (default prefix = "reads_bench"):
    reads_bench_ecoli.png
    reads_bench_human.png
    reads_bench_combined.png
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


TOOLS = ["tuna", "kmc", "fastk", "kfc"]

STYLE = {
    "tuna":  dict(color="#2166ac", label="Tuna"),
    "kmc":   dict(color="#d6604d", label="KMC3"),
    "fastk": dict(color="#4dac26", label="FastK"),
    "kfc":   dict(color="#8856a7", label="KFC"),
}


# ── parsing ──────────────────────────────────────────────────────────────────

def _wall_to_sec(s: str) -> float:
    parts = s.strip().split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def parse_log(path: Path) -> list[dict]:
    hdr_re  = re.compile(r"=== TOOL=(\w+)\s+DS=(\w+)\s+FILE=\S+\s+T=(\d+)")
    wall_re = re.compile(r"Elapsed \(wall clock\) time.*?:\s+([\d:\.]+)")
    rss_re  = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
    ok_re   = re.compile(r"Exit status:\s+0")
    kill_re = re.compile(r"Command terminated by signal")

    text   = path.read_text(errors="replace")
    chunks = re.split(r"(?=^=== TOOL=)", text, flags=re.MULTILINE)

    records = []
    for chunk in chunks:
        m = hdr_re.search(chunk)
        if not m:
            continue
        tool    = m.group(1).lower()
        dataset = m.group(2).lower()
        threads = int(m.group(3))

        wm  = wall_re.search(chunk)
        rm  = rss_re.search(chunk)
        ok  = bool(ok_re.search(chunk)) and not bool(kill_re.search(chunk))

        records.append(dict(
            tool    = tool,
            dataset = dataset,
            threads = threads,
            wall_s  = _wall_to_sec(wm.group(1)) if wm else None,
            rss_kb  = int(rm.group(1))           if rm else None,
            ok      = ok,
        ))
    return records


def resolve_path(p: str) -> Path:
    path = Path(p)
    if path.is_dir():
        candidate = path / "results.log"
        if candidate.exists():
            return candidate
        sys.exit(f"error: no results.log in {path}")
    return path


# ── formatters ────────────────────────────────────────────────────────────────

def _fmt_time(v, _):
    if v <= 0:
        return "0 s"
    if v < 60:
        return f"{v:.0f} s"
    return f"{v/60:.0f} min"


def _fmt_rss(v, _, unit_gb=True):
    if unit_gb:
        return f"{v:.1f} GB"
    return f"{v:.0f} MB"


# ── core bar-chart drawing ────────────────────────────────────────────────────

def draw_time_bars(ax, all_records: list[tuple[str, list]], dataset: str,
                   logscale: bool = False, add_legend: bool = True):
    """Grouped bar chart: tools on x-axis, one bar group per run label."""
    n_runs  = len(all_records)
    n_tools = len(TOOLS)
    width   = 0.7 / max(n_runs, 1)
    offsets = np.linspace(-(n_runs - 1) / 2, (n_runs - 1) / 2, n_runs) * width

    x = np.arange(n_tools)
    any_crash = False

    for ri, (label, records) in enumerate(all_records):
        for ti, tool in enumerate(TOOLS):
            hits = [r for r in records
                    if r["tool"] == tool and r["dataset"] == dataset]
            if not hits:
                continue
            r    = hits[0]
            val  = r["wall_s"]
            ok   = r["ok"]
            if val is None:
                continue

            color   = STYLE[tool]["color"]
            hatch   = "///" if not ok else None
            bar_lbl = label if ti == 0 else None
            ax.bar(x[ti] + offsets[ri], val, width,
                   color=color, alpha=(0.55 if n_runs > 1 and ri > 0 else 0.9),
                   hatch=hatch, edgecolor="white", linewidth=0.6,
                   label=bar_lbl)
            if not ok:
                any_crash = True

    ax.set_xticks(x)
    ax.set_xticklabels([STYLE[t]["label"] for t in TOOLS], fontsize=9)
    ax.set_ylabel("Wall-clock time", fontsize=9)

    if logscale:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    else:
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))

    if any_crash:
        ax.annotate("/// = killed / crashed", xy=(0.97, 0.97),
                    xycoords="axes fraction", ha="right", va="top",
                    fontsize=7, color="grey", style="italic")

    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)

    if add_legend and n_runs > 1:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper right")


def draw_rss_bars(ax, all_records: list[tuple[str, list]], dataset: str,
                  unit_gb: bool = True, add_legend: bool = False):
    n_runs  = len(all_records)
    n_tools = len(TOOLS)
    width   = 0.7 / max(n_runs, 1)
    offsets = np.linspace(-(n_runs - 1) / 2, (n_runs - 1) / 2, n_runs) * width
    unit_div   = 1024**2 if unit_gb else 1024
    unit_label = "GB" if unit_gb else "MB"

    x = np.arange(n_tools)

    for ri, (label, records) in enumerate(all_records):
        for ti, tool in enumerate(TOOLS):
            hits = [r for r in records
                    if r["tool"] == tool and r["dataset"] == dataset]
            if not hits:
                continue
            r   = hits[0]
            val = r["rss_kb"]
            if val is None:
                continue
            color   = STYLE[tool]["color"]
            bar_lbl = label if ti == 0 else None
            ax.bar(x[ti] + offsets[ri], val / unit_div, width,
                   color=color, alpha=(0.55 if n_runs > 1 and ri > 0 else 0.9),
                   edgecolor="white", linewidth=0.6,
                   label=bar_lbl)

    ax.set_xticks(x)
    ax.set_xticklabels([STYLE[t]["label"] for t in TOOLS], fontsize=9)
    ax.set_ylabel(f"Peak RSS ({unit_label})", fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)

    if add_legend and n_runs > 1:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper right")


# ── figures ──────────────────────────────────────────────────────────────────

DATASET_META = {
    "ecoli": dict(title="E. coli reads (k=31, SRR2584863, ~50x)",
                  logtime=False, rss_gb=False),
    "human": dict(title="Human reads (k=31, SRR622461, NA12878 ~8x)",
                  logtime=True, rss_gb=True),
}


def plot_dataset(all_records, dataset, meta, outprefix):
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9, 3.8),
                                      constrained_layout=True)
    fig.suptitle(meta["title"], fontsize=11, fontweight="bold")

    draw_time_bars(ax_t, all_records, dataset,
                   logscale=meta["logtime"], add_legend=True)
    draw_rss_bars(ax_r, all_records, dataset,
                  unit_gb=meta["rss_gb"], add_legend=False)

    out = f"{outprefix}_{dataset}.png"
    fig.savefig(out, dpi=500)
    plt.close(fig)
    print(f"  saved  {out}")


def plot_combined(all_records, outprefix):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7),
                              constrained_layout=True)
    fig.suptitle("k-mer counting benchmark — sequencing reads  (k=31, T=8)",
                 fontsize=11, fontweight="bold")

    datasets   = ["ecoli", "human"]
    row_labels = ["E. coli", "Human"]

    for row, (ds, rl) in enumerate(zip(datasets, row_labels)):
        meta = DATASET_META[ds]

        draw_time_bars(axes[row, 0], all_records, ds,
                       logscale=meta["logtime"],
                       add_legend=(row == 0))
        draw_rss_bars(axes[row, 1], all_records, ds,
                      unit_gb=meta["rss_gb"],
                      add_legend=False)

        time_ylabel = "Wall-clock time (log)" if meta["logtime"] else "Wall-clock time"
        axes[row, 0].set_ylabel(f"{rl}\n{time_ylabel}", fontsize=9)
        axes[row, 1].set_ylabel(
            f"Peak RSS ({'GB' if meta['rss_gb'] else 'MB'})", fontsize=9)

    out = f"{outprefix}_combined.png"
    fig.savefig(out, dpi=500)
    plt.close(fig)
    print(f"  saved  {out}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="results.log files or dirs containing results.log")
    ap.add_argument("--label", nargs="+", metavar="L",
                    help="one label per input (default: parent dir name)")
    ap.add_argument("--out", default="reads_bench",
                    help="output file prefix (default: reads_bench)")
    args = ap.parse_args()

    paths = [resolve_path(p) for p in args.inputs]

    if args.label:
        if len(args.label) != len(paths):
            sys.exit("error: --label count must match number of inputs")
        labels = args.label
    else:
        labels = [p.parent.name for p in paths]

    all_records: list[tuple[str, list]] = []
    for label, path in zip(labels, paths):
        if not path.exists():
            sys.exit(f"error: {path} not found")
        recs = parse_log(path)
        if not recs:
            sys.exit(f"error: no run blocks found in {path}")
        print(f"Parsed {len(recs)} run(s) from {path}  [{label}]")
        all_records.append((label, recs))

    for ds, meta in DATASET_META.items():
        plot_dataset(all_records, ds, meta, args.out)

    plot_combined(all_records, args.out)


if __name__ == "__main__":
    main()
