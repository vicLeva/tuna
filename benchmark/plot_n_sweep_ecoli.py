#!/usr/bin/env python3
"""Plot n_sweep_ecoli.tsv: phase1/phase2/total wall time and RSS vs partition count."""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

tsv = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path(__file__).parent / "results" / "n_sweep_ecoli.tsv"

df = pd.read_csv(tsv, sep="\t")

stats = df.groupby("n").agg(
    p1_mean=("phase1", "mean"), p1_std=("phase1", "std"),
    p2_mean=("phase2", "mean"), p2_std=("phase2", "std"),
    tot_mean=("total",  "mean"), tot_std=("total",  "std"),
    rss_mean=("rss_mb", "mean"), rss_std=("rss_mb", "std"),
).reset_index()

ns = stats["n"].values

C_P1  = "#4e79a7"
C_P2  = "#76b7b2"
C_TOT = "#f28e2b"
C_RSS = "#e15759"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("tuna — E. coli 200 files, k=31  (n sweep)", fontweight="bold")

# ── Wall time ─────────────────────────────────────────────────────────────────
for mean, std, color, label in [
    (stats["p1_mean"],  stats["p1_std"],  C_P1,  "phase1 (partition)"),
    (stats["p2_mean"],  stats["p2_std"],  C_P2,  "phase2 (count+write)"),
    (stats["tot_mean"], stats["tot_std"], C_TOT, "total"),
]:
    ax1.plot(ns, mean, marker="o", color=color, linewidth=2, label=label)
    ax1.fill_between(ns, mean - std, mean + std, alpha=0.15, color=color)
    for n, m in zip(ns, mean):
        ax1.annotate(f"{m:.1f}s", (n, m), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=7, color=color)

ax1.set_xscale("log", base=2)
ax1.set_xticks(ns)
ax1.set_xticklabels([str(n) for n in ns], rotation=45, ha="right", fontsize=8)
ax1.xaxis.set_minor_formatter(plt.NullFormatter())
ax1.set_xlabel("n partitions")
ax1.set_ylabel("wall time (s)")
ax1.set_title("Wall time vs n")
ax1.legend(fontsize=9)
ax1.grid(which="both", alpha=0.25)

# ── RSS ───────────────────────────────────────────────────────────────────────
ax2.plot(ns, stats["rss_mean"], marker="s", color=C_RSS, linewidth=2, label="RSS")
ax2.fill_between(ns, stats["rss_mean"] - stats["rss_std"],
                     stats["rss_mean"] + stats["rss_std"], alpha=0.15, color=C_RSS)
for n, m in zip(ns, stats["rss_mean"]):
    ax2.annotate(f"{m:.0f}", (n, m), textcoords="offset points",
                 xytext=(0, 6), ha="center", fontsize=7, color=C_RSS)

ax2.set_xscale("log", base=2)
ax2.set_xticks(ns)
ax2.set_xticklabels([str(n) for n in ns], rotation=45, ha="right", fontsize=8)
ax2.xaxis.set_minor_formatter(plt.NullFormatter())
ax2.set_xlabel("n partitions")
ax2.set_ylabel("RSS (MB)")
ax2.set_title("Peak memory vs n")
ax2.legend(fontsize=9)
ax2.grid(which="both", alpha=0.25)

out = tsv.parent / "n_sweep_ecoli.png"
plt.tight_layout()
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved: {out}")
plt.show()
