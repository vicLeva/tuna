#!/usr/bin/env python3
"""Compare two bench_datasets runs (e.g. tuna-dev vs tuna-ud).

Usage:
    python plot_datasets_compare.py <csv_a> <label_a> <csv_b> <label_b> [out.png]
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from collections import defaultdict

# ── Args ──────────────────────────────────────────────────────────────────────

if len(sys.argv) < 5:
    print(__doc__)
    sys.exit(1)

csv_a   = Path(sys.argv[1])
label_a = sys.argv[2]
csv_b   = Path(sys.argv[3])
label_b = sys.argv[4]
out_path = Path(sys.argv[5]) if len(sys.argv) > 5 else csv_a.parent / "compare.png"

# ── Load ──────────────────────────────────────────────────────────────────────

def load(path):
    tuna = defaultdict(dict)
    kmc  = defaultdict(dict)
    def _f(row, col):
        v = row.get(col, "")
        return float(v) if v not in ("", "na") else 0.0
    with open(path) as f:
        for row in csv.DictReader(f):
            ds, fidx, tool = row["dataset"], int(row["file_idx"]), row["tool"]
            d = dict(wall=_f(row,"wall_s"), p1=_f(row,"phase1_s"),
                     p2=_f(row,"phase2_s"), rss=_f(row,"rss_mb"))
            if tool == "tuna":
                tuna[ds][fidx] = d
            elif tool == "kmc":
                kmc[ds][fidx] = d
    return tuna, kmc

tuna_a, kmc_a = load(csv_a)
tuna_b, kmc_b = load(csv_b)

DATASETS = [ds for ds in ["ecoli", "salmonella", "gut", "human", "tara"]
            if ds in tuna_a and ds in tuna_b]
ND = len(DATASETS)

# ── Colours ───────────────────────────────────────────────────────────────────

C_A    = "#4e79a7"   # dev  (blue)
C_A2   = "#76b7b2"   # dev  phase2 accent
C_B    = "#9467bd"   # ud   (purple)
C_B2   = "#c5b0d5"   # ud   phase2 accent
C_KMC  = "#f28e2b"   # KMC  (orange)
C_WIN  = "#59a14f"   # speedup > 1 (green)
C_LOSE = "#e15759"   # speedup < 1 (red)

# ── Figure: 3 rows × ND cols ─────────────────────────────────────────────────
# Row 0 — wall time box plots: label_a / label_b / KMC
# Row 1 — phase breakdown: side-by-side stacked bars label_a vs label_b
# Row 2 — per-file relative speedup: label_a_wall / label_b_wall

fig, axes = plt.subplots(3, ND, figsize=(4 * ND, 13))
if ND == 1:
    axes = axes.reshape(3, 1)

fig.suptitle(f"{label_a} vs {label_b} — per-file benchmark  (k=31, m=21, 8 threads)",
             fontsize=13, fontweight="bold")

for col, ds in enumerate(DATASETS):
    common = sorted(set(tuna_a[ds]) & set(tuna_b[ds]) & set(kmc_a[ds]))
    wa  = np.array([tuna_a[ds][i]["wall"] for i in common])
    wb  = np.array([tuna_b[ds][i]["wall"] for i in common])
    wk  = np.array([kmc_a [ds][i]["wall"] for i in common])
    p1a = np.array([tuna_a[ds][i]["p1"]   for i in common])
    p2a = np.array([tuna_a[ds][i]["p2"]   for i in common])
    p1b = np.array([tuna_b[ds][i]["p1"]   for i in common])
    p2b = np.array([tuna_b[ds][i]["p2"]   for i in common])
    n   = len(common)

    # ── Row 0: wall time box plots ─────────────────────────────────────────
    ax = axes[0, col]
    bp = ax.boxplot([wa, wb, wk], patch_artist=True,
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(linewidth=1.2),
                    flierprops=dict(marker=".", markersize=4, alpha=0.5),
                    widths=0.5)
    for box, c in zip(bp["boxes"], [C_A, C_B, C_KMC]):
        box.set_facecolor(c)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([label_a, label_b, "KMC"], fontsize=8)
    ax.set_title(f"{ds}  (n={n})", fontweight="bold")
    ax.set_ylabel("Wall time (s)" if col == 0 else "")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    for pos, vals, c in [(1, wa, C_A), (2, wb, C_B), (3, wk, C_KMC)]:
        med = np.median(vals)
        ax.text(pos, ax.get_ylim()[1] * 0.97, f"{med:.2f}s",
                ha="center", va="top", fontsize=7, color=c, fontweight="bold")

    # ── Row 1: phase breakdown — side-by-side stacked bars ────────────────
    ax = axes[1, col]
    order = np.argsort(p1a + p2a)
    x = np.arange(n)
    w = 0.42
    # label_a stack
    ax.bar(x - w/2, p1a[order], width=w, color=C_A,  label=f"{label_a} p1" if col == 0 else "")
    ax.bar(x - w/2, p2a[order], width=w, color=C_A2, bottom=p1a[order],
           label=f"{label_a} p2" if col == 0 else "")
    # label_b stack (sorted by same order for comparability)
    ax.bar(x + w/2, p1b[order], width=w, color=C_B,  label=f"{label_b} p1" if col == 0 else "")
    ax.bar(x + w/2, p2b[order], width=w, color=C_B2, bottom=p1b[order],
           label=f"{label_b} p2" if col == 0 else "")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("Files (sorted by total)")
    ax.set_ylabel("Wall time (s)" if col == 0 else "")
    ax.set_title("Phase breakdown")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    if col == 0:
        ax.legend(fontsize=7, loc="upper left", ncol=2)

    # ── Row 2: relative speedup label_a / label_b per file ────────────────
    ax = axes[2, col]
    ratio = wa / wb   # > 1 means label_b faster
    order2 = np.argsort(ratio)
    colors = [C_WIN if r >= 1 else C_LOSE for r in ratio[order2]]
    ax.bar(np.arange(n), ratio[order2], color=colors, width=1.0)
    ax.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
    med = np.median(ratio)
    ax.axhline(med, color=C_WIN if med >= 1 else C_LOSE, linewidth=1.5, linestyle=":",
               label=f"median {med:.2f}×")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("Files (sorted by ratio)")
    ax.set_ylabel(f"{label_a} / {label_b}  (>1 = {label_b} faster)" if col == 0 else "")
    ax.set_title(f"Speedup  med={med:.2f}×")
    ax.legend(fontsize=8, loc="lower right")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
