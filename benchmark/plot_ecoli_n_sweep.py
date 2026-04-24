#!/usr/bin/env python3
"""
Plot E. coli N-sweep benchmark: tuna vs KMC wall time across thread counts.

CSV columns: tool, n_files, threads, rep, wall_s, rss_mb, phase1_s, phase2_s, unique_kmers
Usage: python plot_ecoli_n_sweep.py <bench_n_sweep.csv> [output.png]
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    sorted(Path(__file__).parent.glob("results/ecoli_nsweep/*/bench_n_sweep.csv"))[-1]
out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else \
    Path(csv_path).parent / "ecoli_n_sweep.png"

df = pd.read_csv(csv_path)
for c in ["wall_s", "rss_mb", "phase1_s", "phase2_s", "unique_kmers"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

thread_counts = sorted(df["threads"].unique())
med = df.groupby(["tool", "n_files", "threads"])[["wall_s", "rss_mb", "phase1_s", "phase2_s"]].median().reset_index()

C_TUNA = "#2196F3"
C_KMC  = "#F44336"

n_cols = len(thread_counts)
fig, axes = plt.subplots(2, n_cols, figsize=(7 * n_cols, 10), squeeze=False)
fig.suptitle(f"tuna vs KMC — E. coli N-sweep  (k=31, m=21)", fontsize=14, fontweight="bold")

for col, t in enumerate(thread_counts):
    sub = med[med["threads"] == t]
    ns = sorted(sub["n_files"].unique())

    # ── Top: wall time ────────────────────────────────────────────────────────
    ax = axes[0][col]
    for tool, color in [("tuna", C_TUNA), ("kmc", C_KMC)]:
        d = sub[sub["tool"] == tool].sort_values("n_files")
        if d.empty: continue
        ax.plot(d["n_files"], d["wall_s"], "o-", color=color, lw=2.5, ms=7,
                label=tool.upper())
        if tool == "tuna":
            ax.plot(d["n_files"], d["phase1_s"], "s--", color=color, lw=1.2, ms=4,
                    alpha=0.55, label="tuna phase1")
            ax.plot(d["n_files"], d["phase2_s"], "^--", color=color, lw=1.2, ms=4,
                    alpha=0.55, label="tuna phase2")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(ns); ax.set_xticklabels([str(n) for n in ns], rotation=45)
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel("N genomes"); ax.set_ylabel("Wall time (s)")
    ax.set_title(f"Wall time — t={t} thread{'s' if t > 1 else ''}")
    ax.legend(fontsize=8); ax.grid(which="both", alpha=0.25)

    # ── Bottom: speedup KMC/tuna ───────────────────────────────────────────────
    ax = axes[1][col]
    tuna_d = sub[sub["tool"] == "tuna"].sort_values("n_files").set_index("n_files")
    kmc_d  = sub[sub["tool"] == "kmc" ].sort_values("n_files").set_index("n_files")
    common = sorted(set(tuna_d.index) & set(kmc_d.index))
    if common:
        speedup = kmc_d.loc[common, "wall_s"].values / tuna_d.loc[common, "wall_s"].values
        ax.axhline(1.0, color="gray", ls="--", lw=1.2)
        ax.plot(common, speedup, "D-", color="#9C27B0", lw=2, ms=7)
        for n, s in zip(common, speedup):
            ax.annotate(f"{s:.2f}×", (n, s), textcoords="offset points",
                        xytext=(0, 6 if s >= 1 else -12), ha="center", fontsize=8,
                        color=C_TUNA if s > 1 else C_KMC, fontweight="bold")
        ax.fill_between(common, 1, speedup,
                         where=speedup >= 1, alpha=0.1, color=C_TUNA, label="tuna faster")
        ax.fill_between(common, speedup, 1,
                         where=speedup < 1,  alpha=0.1, color=C_KMC,  label="KMC faster")

    ax.set_xscale("log")
    ax.set_xticks(common if common else ns)
    ax.set_xticklabels([str(n) for n in (common if common else ns)], rotation=45)
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel("N genomes"); ax.set_ylabel("KMC wall / tuna wall")
    ax.set_title(f"Speedup (>1 = tuna faster) — t={t}")
    ax.legend(fontsize=8); ax.grid(which="both", alpha=0.25)

plt.tight_layout()
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
