#!/usr/bin/env python3
"""Plot per-file tuna vs KMC benchmark results from bench_all_datasets.sh.

Usage:
    python plot_all_datasets.py <results_dir> [out.png]

<results_dir> is the directory produced by bench_all_datasets.sh (e.g. after
scp'ing it back from the cluster) — must contain bench_datasets.csv.
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from collections import defaultdict

# ── Load ──────────────────────────────────────────────────────────────────────

if len(sys.argv) < 2:
    sys.exit("usage: plot_all_datasets.py <results_dir> [out.png]")

results_dir = Path(sys.argv[1])
csv_path = results_dir / "bench_datasets.csv"
if not csv_path.is_file():
    sys.exit(f"error: not found: {csv_path}")

out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else results_dir / "datasets.png"


def _f(v):
    # phase1_s/phase2_s carry a trailing "s" from tuna/KMC's own stderr
    # formatting (e.g. "0.053951s") that the bench script's se_val helper
    # doesn't always strip (only does so when preceded by a space).
    v = v.rstrip("s") if v not in ("", "na") else v
    return float(v) if v not in ("", "na") else None


# tuna[ds][fidx] = {wall, p1, p2, rss, unique}   kmc[ds][fidx] = {...}
tuna = defaultdict(dict)
kmc = defaultdict(dict)

with open(csv_path) as f:
    for row in csv.DictReader(f):
        ds = row["dataset"]
        fidx = int(row["file_idx"])
        entry = dict(
            wall=_f(row["wall_s"]),
            p1=_f(row["phase1_s"]),
            p2=_f(row["phase2_s"]),
            rss=_f(row["rss_mb"]),
            unique=_f(row["unique_kmers"]),
        )
        if row["tool"] == "tuna":
            tuna[ds][fidx] = entry
        elif row["tool"] == "kmc":
            kmc[ds][fidx] = entry

DATASETS = [ds for ds in ["ecoli", "salmonella", "gut", "human", "tara"] if ds in tuna]
ND = len(DATASETS)

# ── Correctness check (console only, never plotted) ──────────────────────────

mismatches = []
for ds in DATASETS:
    common = sorted(set(tuna[ds]) & set(kmc[ds]))
    for fidx in common:
        tu, ku = tuna[ds][fidx]["unique"], kmc[ds][fidx]["unique"]
        if tu is not None and ku is not None and tu != ku:
            mismatches.append((ds, fidx, tu, ku))

if mismatches:
    print(f"[WARN] {len(mismatches)} tuna/KMC unique-kmer mismatch(es):")
    for ds, fidx, tu, ku in mismatches:
        print(f"    {ds} file_idx={fidx}: tuna={tu:.0f}  kmc={ku:.0f}")
else:
    print("[OK] tuna and KMC agree on unique k-mer counts for all common files")

# ── Colours ───────────────────────────────────────────────────────────────────

C_TUNA = "#4e79a7"
C_KMC = "#f28e2b"
C_SPD = "#59a14f"
C_P1 = "#a6c8e8"   # phase1 (lighter shade, tool-agnostic)
C_P2 = "#2c5f8a"   # phase2 (darker shade, tool-agnostic)

# ── Figure layout: 3 rows × ND cols ──────────────────────────────────────────
# Row 0 — wall time distributions (box plots, tuna vs KMC)
# Row 1 — phase1/phase2 breakdown (stacked bars of per-tool medians)
# Row 2 — speedup = kmc_wall / tuna_wall per file

fig, axes = plt.subplots(3, ND, figsize=(4 * ND, 13))
if ND == 1:
    axes = axes.reshape(3, 1)

fig.suptitle("tuna vs KMC3 — per-file benchmark  (k=31, m=21, count-only)",
             fontsize=13, fontweight="bold")

for col, ds in enumerate(DATASETS):
    common = sorted(set(tuna[ds]) & set(kmc[ds]))
    t_wall = np.array([tuna[ds][i]["wall"] for i in common])
    k_wall = np.array([kmc[ds][i]["wall"] for i in common])
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
    ax.set_xticklabels(["tuna", "KMC3"])
    ax.set_title(f"{ds}  (n={n})", fontweight="bold")
    ax.set_ylabel("Wall time (s)" if col == 0 else "")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    for pos, vals, c in [(1, t_wall, C_TUNA), (2, k_wall, C_KMC)]:
        med = np.median(vals)
        ax.text(pos, ax.get_ylim()[1] * 0.97, f"{med:.2f}s",
                ha="center", va="top", fontsize=8, color=c, fontweight="bold")

    # ── Row 1: phase1/phase2 breakdown (median, stacked) ───────────────────
    ax = axes[1, col]
    t_p1 = np.array([tuna[ds][i]["p1"] for i in common if tuna[ds][i]["p1"] is not None])
    t_p2 = np.array([tuna[ds][i]["p2"] for i in common if tuna[ds][i]["p2"] is not None])
    k_p1 = np.array([kmc[ds][i]["p1"] for i in common if kmc[ds][i]["p1"] is not None])
    k_p2 = np.array([kmc[ds][i]["p2"] for i in common if kmc[ds][i]["p2"] is not None])

    med_t_p1 = np.median(t_p1) if len(t_p1) else 0
    med_t_p2 = np.median(t_p2) if len(t_p2) else 0
    med_k_p1 = np.median(k_p1) if len(k_p1) else 0
    med_k_p2 = np.median(k_p2) if len(k_p2) else 0

    xs = [0, 1]
    ax.bar(xs, [med_t_p1, med_k_p1], color=C_P1, width=0.5, label="phase1")
    ax.bar(xs, [med_t_p2, med_k_p2], bottom=[med_t_p1, med_k_p1],
           color=C_P2, width=0.5, label="phase2")
    ax.set_xticks(xs)
    ax.set_xticklabels(["tuna", "KMC3"])
    ax.set_ylabel("Median phase time (s)" if col == 0 else "")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)
    if col == 0:
        ax.legend(fontsize=8, loc="upper right")

    # ── Row 2: speedup distribution ────────────────────────────────────────
    ax = axes[2, col]
    order2 = np.argsort(speedup)
    colors = [C_SPD if s >= 1 else C_KMC for s in speedup[order2]]
    ax.bar(np.arange(n), speedup[order2], color=colors, width=1.0)
    ax.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("Files (sorted by speedup)")
    ax.set_ylabel("KMC3 / tuna  (>1 = tuna faster)" if col == 0 else "")
    ax.set_title(f"Speedup  med={np.median(speedup):.2f}×")
    med = np.median(speedup)
    ax.axhline(med, color=C_SPD, linewidth=1.5, linestyle=":",
               label=f"median {med:.2f}×")
    ax.legend(fontsize=8, loc="lower right")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(out_path, dpi=500, bbox_inches="tight")
print(f"Saved: {out_path}")
