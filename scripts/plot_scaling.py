#!/usr/bin/env python3
"""Plot n-files scaling sweep from bench_scaling.sh output.

Usage:
    python scripts/plot_scaling.py [bench_scaling.csv] [outprefix]

Defaults:
    csv     — most recent benchmark/results/bench_scaling_*/bench_scaling.csv
    outprefix — same directory as csv, stem = "scaling"
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path
from collections import defaultdict

# ── Input paths ───────────────────────────────────────────────────────────────

if len(sys.argv) > 1:
    csv_path = Path(sys.argv[1])
else:
    candidates = sorted(
        Path(__file__).parent.parent.glob(
            "benchmark/results/bench_scaling_*/bench_scaling.csv"
        )
    )
    if not candidates:
        raise FileNotFoundError("No bench_scaling.csv found under benchmark/results/")
    csv_path = candidates[-1]

outprefix = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_path.parent / "scaling"

print(f"CSV:    {csv_path}")
print(f"Output: {outprefix}_*.png")

# ── Style (matches plot_all_tools.py) ─────────────────────────────────────────

STYLE = {
    "tuna":  dict(color="#2166ac", marker="o",  label="tuna"),
    "kmc":   dict(color="#d6604d", marker="s",  label="KMC3"),
    "fastk": dict(color="#4dac26", marker="^",  label="FastK"),
}

# ── Load CSV ──────────────────────────────────────────────────────────────────
# data[dataset][tool][n_files] = {wall_s, rss_mb, phase1_s, phase2_s}

data = defaultdict(lambda: defaultdict(dict))

with open(csv_path) as f:
    for row in csv.DictReader(f):
        ds   = row["dataset"]
        tool = row["tool"]
        n    = int(row["n_files"])
        def _f(k):
            v = row.get(k, "")
            return float(v) if v not in ("", "na") else None
        data[ds][tool][n] = dict(
            wall=_f("wall_s"),
            rss =_f("rss_mb"),
            p1  =_f("phase1_s"),
            p2  =_f("phase2_s"),
        )

DATASETS = [ds for ds in ["ecoli", "human"] if ds in data]
TOOLS    = [t  for t  in ["tuna", "kmc", "fastk"] if any(t in data[ds] for ds in DATASETS)]

DS_LABELS = {"ecoli": "E. coli (assemblies)", "human": "Human (assemblies)"}

# ── Helpers ───────────────────────────────────────────────────────────────────

def xs_ys(ds, tool, metric="wall"):
    pts = data[ds].get(tool, {})
    pairs = sorted((n, v[metric]) for n, v in pts.items() if v.get(metric) is not None)
    if not pairs:
        return [], []
    xs, ys = zip(*pairs)
    return list(xs), list(ys)

# ── Figure: wall time ─────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, len(DATASETS), figsize=(6 * len(DATASETS), 4.5), sharey=False)
if len(DATASETS) == 1:
    axes = [axes]

for ax, ds in zip(axes, DATASETS):
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(ds, tool, "wall")
        if not xs:
            continue
        ax.plot(xs, ys, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=1.8, markersize=5, zorder=3)

    ax.set_xscale("log")
    ax.set_xlabel("Number of input files", fontsize=10)
    ax.set_ylabel("Wall-clock time (s)", fontsize=10)
    ax.set_title(DS_LABELS.get(ds, ds), fontsize=11, fontweight="bold")
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.grid(axis="both", linestyle=":", alpha=0.4)
    ax.spines[["top", "right"]].set_visible(False)

# shared legend
handles, labels = axes[0].get_legend_handles_labels()
seen = set()
h2, l2 = [], []
for h, l in zip(handles, labels):
    if l not in seen:
        h2.append(h); l2.append(l); seen.add(l)

fig.legend(h2, l2, loc="upper left", bbox_to_anchor=(0.01, 0.99),
           fontsize=9, framealpha=0.85)
fig.suptitle("Scaling vs number of input files  (k=31, m=21, 8 threads)",
             fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(f"{outprefix}_wall.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {outprefix}_wall.png")

# ── Figure: phase breakdown (tuna only) ───────────────────────────────────────

if "tuna" in TOOLS:
    fig, axes = plt.subplots(1, len(DATASETS), figsize=(6 * len(DATASETS), 4.5), sharey=False)
    if len(DATASETS) == 1:
        axes = [axes]

    for ax, ds in zip(axes, DATASETS):
        xs, p1s = xs_ys(ds, "tuna", "p1")
        _,  p2s = xs_ys(ds, "tuna", "p2")
        if xs:
            ax.plot(xs, p1s, color="#4393c3", marker="o", linewidth=1.6,
                    markersize=5, label="phase 1 (partition)")
            ax.plot(xs, p2s, color="#d6604d", marker="s", linewidth=1.6,
                    markersize=5, label="phase 2 (count+write)")
            ax.plot(xs, [a + b for a, b in zip(p1s, p2s)],
                    color="#2166ac", marker="o", linewidth=1.6, markersize=5,
                    linestyle="--", label="total", zorder=2)

        ax.set_xscale("log")
        ax.set_xlabel("Number of input files", fontsize=10)
        ax.set_ylabel("Time (s)", fontsize=10)
        ax.set_title(f"tuna — {DS_LABELS.get(ds, ds)}", fontsize=11, fontweight="bold")
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.grid(axis="both", linestyle=":", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8.5)

    fig.suptitle("tuna phase breakdown vs number of input files  (k=31, m=21, 8 threads)",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(f"{outprefix}_phases.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outprefix}_phases.png")

# ── Figure: peak RSS ──────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, len(DATASETS), figsize=(6 * len(DATASETS), 4.5), sharey=False)
if len(DATASETS) == 1:
    axes = [axes]

for ax, ds in zip(axes, DATASETS):
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = xs_ys(ds, tool, "rss")
        if not xs:
            continue
        ys_gb = [v / 1024 for v in ys]
        ax.plot(xs, ys_gb, color=st["color"], marker=st["marker"],
                label=st["label"], linewidth=1.8, markersize=5, zorder=3)

    ax.set_xscale("log")
    ax.set_xlabel("Number of input files", fontsize=10)
    ax.set_ylabel("Peak RSS (GB)", fontsize=10)
    ax.set_title(DS_LABELS.get(ds, ds), fontsize=11, fontweight="bold")
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.grid(axis="both", linestyle=":", alpha=0.4)
    ax.spines[["top", "right"]].set_visible(False)

handles, labels = axes[0].get_legend_handles_labels()
seen = set()
h2, l2 = [], []
for h, l in zip(handles, labels):
    if l not in seen:
        h2.append(h); l2.append(l); seen.add(l)

fig.legend(h2, l2, loc="upper left", bbox_to_anchor=(0.01, 0.99),
           fontsize=9, framealpha=0.85)
fig.suptitle("Peak RSS vs number of input files  (k=31, m=21, 8 threads)",
             fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(f"{outprefix}_rss.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {outprefix}_rss.png")
