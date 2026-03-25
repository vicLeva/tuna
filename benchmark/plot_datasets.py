#!/usr/bin/env python3
"""Plot per-file tuna vs KMC benchmark results from bench_datasets.sh.

Usage:
    python plot_datasets.py [bench.csv] [out.png]
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from collections import defaultdict

# ── Load ──────────────────────────────────────────────────────────────────────

csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    sorted(Path(__file__).parent.glob("../benchmark/results/bench_datasets_*/bench.csv"))[-1]
out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else \
    Path(csv_path).parent / "datasets.png"

# tuna[ds][fidx] = {wall, p1, p2, rss}   kmc[ds][fidx] = {wall, rss}
tuna = defaultdict(dict)
kmc  = defaultdict(dict)

with open(csv_path) as f:
    for row in csv.DictReader(f):
        ds   = row["dataset"]
        fidx = int(row["file_idx"])
        tool = row["tool"]
        if tool == "tuna":
            tuna[ds][fidx] = dict(
                wall  = float(row["wall_s"]),
                p1    = float(row["phase1_s"]) if row["phase1_s"] not in ("", "na") else 0,
                p2    = float(row["phase2_s"]) if row["phase2_s"] not in ("", "na") else 0,
                rss   = float(row["rss_mb"]),
            )
        elif tool == "kmc":
            kmc[ds][fidx] = dict(
                wall  = float(row["wall_s"]),
                p1    = float(row["phase1_s"]) if row["phase1_s"] not in ("", "na") else 0,
                p2    = float(row["phase2_s"]) if row["phase2_s"] not in ("", "na") else 0,
                rss   = float(row["rss_mb"]),
            )

DATASETS = [ds for ds in ["ecoli", "human", "salmonella", "gut", "tara"] if ds in tuna]
ND = len(DATASETS)

# ── Colours ───────────────────────────────────────────────────────────────────

C_TUNA = "#4e79a7"
C_KMC  = "#f28e2b"
C_P1   = "#4e79a7"
C_P2   = "#76b7b2"
C_SPD  = "#59a14f"

# ── Figure layout: 3 rows × ND cols ──────────────────────────────────────────
# Row 0 — wall time distributions (box plots, tuna vs KMC)
# Row 1 — tuna phase breakdown (p1 / p2 per file, sorted by total)
# Row 2 — speedup = kmc_wall / tuna_wall per file

fig, axes = plt.subplots(3, ND, figsize=(4 * ND, 13))
if ND == 1:
    axes = axes.reshape(3, 1)

fig.suptitle("tuna vs KMC — per-file benchmark  (k=31, m=21, 8 threads)",
             fontsize=13, fontweight="bold")

for col, ds in enumerate(DATASETS):
    common = sorted(set(tuna[ds]) & set(kmc[ds]))
    t_wall = np.array([tuna[ds][i]["wall"] for i in common])
    k_wall = np.array([kmc [ds][i]["wall"] for i in common])
    t_p1   = np.array([tuna[ds][i]["p1"]   for i in common])
    t_p2   = np.array([tuna[ds][i]["p2"]   for i in common])
    speedup = k_wall / t_wall

    n = len(common)

    # ── Row 0: wall time box plots ─────────────────────────────────────────
    ax = axes[0, col]
    bp = ax.boxplot([t_wall, k_wall], patch_artist=True,
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(linewidth=1.2),
                    flierprops=dict(marker=".", markersize=4, alpha=0.5),
                    widths=0.5)
    bp["boxes"][0].set_facecolor(C_TUNA)
    bp["boxes"][1].set_facecolor(C_KMC)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["tuna", "KMC"])
    ax.set_title(f"{ds}  (n={n})", fontweight="bold")
    ax.set_ylabel("Wall time (s)" if col == 0 else "")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    # annotate medians
    for pos, vals, c in [(1, t_wall, C_TUNA), (2, k_wall, C_KMC)]:
        med = np.median(vals)
        ax.text(pos, ax.get_ylim()[1] * 0.97, f"{med:.2f}s",
                ha="center", va="top", fontsize=8, color=c, fontweight="bold")

    # ── Row 1: tuna phase breakdown sorted by total ────────────────────────
    ax = axes[1, col]
    order = np.argsort(t_p1 + t_p2)
    x = np.arange(n)
    ax.bar(x, t_p1[order], color=C_P1, label="phase1" if col == 0 else "", width=1.0)
    ax.bar(x, t_p2[order], bottom=t_p1[order], color=C_P2,
           label="phase2" if col == 0 else "", width=1.0)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("Files (sorted by total)")
    ax.set_ylabel("tuna wall time (s)" if col == 0 else "")
    ax.set_title("tuna phase breakdown")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    if col == 0:
        ax.legend(fontsize=8, loc="upper left")

    # ── Row 2: speedup distribution ────────────────────────────────────────
    ax = axes[2, col]
    order2 = np.argsort(speedup)
    colors = [C_SPD if s >= 1 else C_KMC for s in speedup[order2]]
    ax.bar(np.arange(n), speedup[order2], color=colors, width=1.0)
    ax.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("Files (sorted by speedup)")
    ax.set_ylabel("KMC / tuna  (>1 = tuna faster)" if col == 0 else "")
    ax.set_title(f"Speedup  med={np.median(speedup):.2f}×")
    med = np.median(speedup)
    ax.axhline(med, color=C_SPD, linewidth=1.5, linestyle=":",
               label=f"median {med:.2f}×")
    ax.legend(fontsize=8, loc="lower right")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
plt.show()
