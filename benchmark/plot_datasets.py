#!/usr/bin/env python3
"""Plot bench_datasets.sh results: one figure per dataset + one summary figure.

Usage:
    python plot_datasets.py <bench.csv>

Outputs (next to bench.csv):
    <dataset>.png   — per-dataset figure (wall time, phase breakdown, RSS, variability)
    summary.png     — cross-dataset comparison (tuna best-m vs KMC)

CSV schema (bench_datasets.sh):
    dataset, file_idx, filename, tool, m, wall_s, rss_mb,
    unique_kmers, phase1_s, phase2_s

  tuna: m in {17,19,21,23}; phase1=partition, phase2=count+write
  kmc:  m="na";             phase1=kmc_count, phase2=kmc_dump, wall=sum
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from collections import defaultdict

# ── Load ──────────────────────────────────────────────────────────────────────

csv_path = Path(sys.argv[1])
out_dir  = csv_path.parent

rows = defaultdict(list)        # (dataset, tool, m) → [dicts]
datasets_ordered = []
ds_nfiles = defaultdict(set)   # dataset → set of file_idx

with open(csv_path) as f:
    for row in csv.DictReader(f):
        ds   = row["dataset"]
        tool = row["tool"]
        m    = row["m"]
        try:
            wall = float(row["wall_s"])
            rss  = float(row["rss_mb"])
            p1   = float(row["phase1_s"])
            p2   = float(row["phase2_s"])
        except (ValueError, KeyError):
            continue
        rows[(ds, tool, m)].append({"wall": wall, "rss": rss, "p1": p1, "p2": p2})
        ds_nfiles[ds].add(row["file_idx"])
        if ds not in datasets_ordered:
            datasets_ordered.append(ds)

m_values = sorted({int(m) for (_, tool, m) in rows if tool == "tuna"})

# ── Helpers ───────────────────────────────────────────────────────────────────

def vals(ds, tool, m, field):
    return [r[field] for r in rows.get((ds, tool, str(m)), [])]

def med(ds, tool, m, field):
    v = vals(ds, tool, m, field)
    return float(np.median(v)) if v else np.nan

def best_m(ds):
    best, bw = None, np.inf
    for m in m_values:
        w = med(ds, "tuna", m, "wall")
        if not np.isnan(w) and w < bw:
            bw, best = w, m
    return best

# ── Colours ───────────────────────────────────────────────────────────────────

C_KMC    = "#e15759"
C_P1_T   = "#4e79a7"
C_P2_T   = "#f28e2b"
C_P1_K   = "#59a14f"
C_P2_K   = "#b07aa1"
C_RSS_T  = "#4e79a7"
C_RSS_K  = "#e15759"

M_COLORS = {17: "#4e79a7", 19: "#f28e2b", 21: "#59a14f", 23: "#b07aa1"}

# x-axis groups: m=17, m=19, m=21, m=23, KMC
GROUPS   = [str(m) for m in m_values] + ["kmc"]
GLABELS  = [f"m={m}" for m in m_values] + ["KMC"]
x        = np.arange(len(GROUPS))
w        = 0.55

# ── Per-dataset figures ───────────────────────────────────────────────────────

for ds in datasets_ordered:

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    n_files = len(ds_nfiles[ds])
    fig.suptitle(f"{ds}   (k=31, {n_files} file(s), median + per-file spread)",
                 fontsize=12, fontweight="bold")

    # ── Panel 1: Wall time stacked (p1 + p2) ─────────────────────────────────

    ax = axes[0]
    for i, grp in enumerate(GROUPS):
        if grp == "kmc":
            p1v = med(ds, "kmc", "na", "p1")
            p2v = med(ds, "kmc", "na", "p2")
            ax.bar(i, p1v, width=w, color=C_P1_K, zorder=3,
                   label="kmc count"  if i == len(GROUPS)-1 else None)
            ax.bar(i, p2v, width=w, bottom=p1v, color=C_P2_K, zorder=3,
                   label="kmc dump"   if i == len(GROUPS)-1 else None)
            total = p1v + p2v
        else:
            m = int(grp)
            p1v = med(ds, "tuna", m, "p1")
            p2v = med(ds, "tuna", m, "p2")
            ax.bar(i, p1v, width=w, color=C_P1_T, zorder=3,
                   label="tuna p1 (partition)"   if i == 0 else None)
            ax.bar(i, p2v, width=w, bottom=p1v, color=C_P2_T, zorder=3,
                   label="tuna p2 (count+write)" if i == 0 else None)
            total = p1v + p2v
        if not np.isnan(total):
            ax.text(i, total * 1.06, f"{total:.2f}s", ha="center",
                    va="bottom", fontsize=8)

    ax.set_xticks(x); ax.set_xticklabels(GLABELS, fontsize=9)
    ax.set_ylabel("Median wall time (s)")
    ax.set_title("Wall time (phase breakdown)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3, zorder=0)

    # ── Panel 2: Per-file wall time box plot ──────────────────────────────────

    ax = axes[1]
    box_data, box_pos, box_colors = [], [], []
    for i, grp in enumerate(GROUPS):
        if grp == "kmc":
            v = vals(ds, "kmc", "na", "wall")
        else:
            v = vals(ds, "tuna", int(grp), "wall")
        if v:
            box_data.append(v)
            box_pos.append(i)
            box_colors.append(C_KMC if grp == "kmc" else M_COLORS[int(grp)])

    if box_data:
        bp = ax.boxplot(box_data, positions=box_pos, widths=0.45,
                        patch_artist=True, showfliers=True,
                        medianprops=dict(color="black", linewidth=1.8))
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)

    ax.set_xticks(x); ax.set_xticklabels(GLABELS, fontsize=9)
    ax.set_ylabel("Wall time per file (s)")
    ax.set_title("Per-file variability")
    ax.grid(axis="y", alpha=0.3, zorder=0)

    # ── Panel 3: Peak RSS ─────────────────────────────────────────────────────

    ax = axes[2]
    for i, grp in enumerate(GROUPS):
        if grp == "kmc":
            rv = med(ds, "kmc", "na", "rss")
            color = C_RSS_K
        else:
            rv = med(ds, "tuna", int(grp), "rss")
            color = M_COLORS[int(grp)]
        ax.bar(i, rv, width=w, color=color, zorder=3)
        if not np.isnan(rv):
            ax.text(i, rv * 1.06, f"{rv:.0f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x); ax.set_xticklabels(GLABELS, fontsize=9)
    ax.set_ylabel("Median peak RSS (MB)")
    ax.set_title("Peak memory")
    ax.grid(axis="y", alpha=0.3, zorder=0)

    plt.tight_layout()
    out = out_dir / f"{ds}.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close()

# ── Summary figure: tuna best-m vs KMC across all datasets ───────────────────

DS    = datasets_ordered
n_ds  = len(DS)
xs    = np.arange(n_ds)
wb    = 0.35

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle("Summary: tuna (best m) vs KMC — all datasets, k=31, median across files",
             fontsize=12, fontweight="bold")

# Panel 1: wall time
ax = axes[0]
tuna_w, kmc_w, best_ms = [], [], []
for ds in DS:
    bm = best_m(ds)
    tuna_w.append(med(ds, "tuna", bm, "wall") if bm else np.nan)
    kmc_w.append(med(ds, "kmc",  "na", "wall"))
    best_ms.append(bm)
ax.bar(xs - wb/2, tuna_w, width=wb, color=C_P1_T, label="tuna (best m)", zorder=3)
ax.bar(xs + wb/2, kmc_w,  width=wb, color=C_KMC,  label="kmc",           zorder=3)
for xi, tv, kv, bm in zip(xs, tuna_w, kmc_w, best_ms):
    if not np.isnan(tv): ax.text(xi - wb/2, tv * 1.08, f"{tv:.1f}s\nm={bm}",
                                  ha="center", va="bottom", fontsize=8)
    if not np.isnan(kv): ax.text(xi + wb/2, kv * 1.08, f"{kv:.1f}s",
                                  ha="center", va="bottom", fontsize=8)
ax.set_xticks(xs); ax.set_xticklabels(DS, fontsize=9)
ax.set_yscale("log"); ax.set_ylabel("Median wall time (s)  [log]")
ax.set_title("Wall time"); ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3, which="both", zorder=0)
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

# Panel 2: speedup tuna / kmc
ax = axes[1]
speedups = [kv / tv if (not np.isnan(tv) and not np.isnan(kv) and tv > 0) else np.nan
            for tv, kv in zip(tuna_w, kmc_w)]
colors_su = [("#59a14f" if s >= 1 else "#e15759") if not np.isnan(s) else "grey"
             for s in speedups]
ax.bar(xs, speedups, width=0.55, color=colors_su, zorder=3)
ax.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
for xi, s in zip(xs, speedups):
    if not np.isnan(s):
        ax.text(xi, s * 1.04 if s >= 1 else s * 0.93,
                f"{s:.2f}×", ha="center", va="bottom" if s >= 1 else "top", fontsize=9)
ax.set_xticks(xs); ax.set_xticklabels(DS, fontsize=9)
ax.set_ylabel("KMC wall / tuna wall  (>1 = tuna faster)")
ax.set_title("Speedup (tuna best-m vs KMC)")
ax.grid(axis="y", alpha=0.3, zorder=0)

# Panel 3: RSS
ax = axes[2]
tuna_r, kmc_r = [], []
for ds in DS:
    bm = best_m(ds)
    tuna_r.append(med(ds, "tuna", bm, "rss") if bm else np.nan)
    kmc_r.append(med(ds, "kmc", "na", "rss"))
ax.bar(xs - wb/2, tuna_r, width=wb, color=C_RSS_T, label="tuna (best m)", zorder=3)
ax.bar(xs + wb/2, kmc_r,  width=wb, color=C_RSS_K,  label="kmc",           zorder=3)
for xi, tv, kv in zip(xs, tuna_r, kmc_r):
    if not np.isnan(tv): ax.text(xi - wb/2, tv * 1.08, f"{tv:.0f}",
                                  ha="center", va="bottom", fontsize=8)
    if not np.isnan(kv): ax.text(xi + wb/2, kv * 1.08, f"{kv:.0f}",
                                  ha="center", va="bottom", fontsize=8)
ax.set_xticks(xs); ax.set_xticklabels(DS, fontsize=9)
ax.set_yscale("log"); ax.set_ylabel("Median peak RSS (MB)  [log]")
ax.set_title("Peak memory"); ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3, which="both", zorder=0)
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

plt.tight_layout()
out = out_dir / "summary.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
plt.close()
