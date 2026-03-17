#!/usr/bin/env python3
"""Plot m-sweep benchmark results (phase timing, overflow, resize analysis)."""

import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# ── Load data ─────────────────────────────────────────────────────────────────

csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    sorted(Path.home().glob("tuna_bench_tmp/bench_m_human_*/bench_m.csv"))[-1]
title = sys.argv[2] if len(sys.argv) > 2 else \
    "Impact of minimizer length m  (k=31, n=4096, t=4, 1 human genome gz)"
kmc_wall = float(sys.argv[3]) if len(sys.argv) > 3 else None  # optional KMC reference time (s)

m_vals, p1, p2 = [], [], []
ov_pct, ov_kmers = [], []
n_resizes, resize_s, n_ov_res, n_load_res = [], [], [], []

with open(csv_path) as f:
    for row in csv.DictReader(f):
        m_vals.append(int(row["m"]))
        p1.append(float(row["phase1_s"]))
        p2.append(float(row["phase2_s"]))
        ov_pct.append(float(row["overflow_pct"]))
        ov_kmers.append(int(row["overflow_kmers"]) if row["overflow_kmers"] else 0)
        n_resizes.append(int(row["n_resizes"])   if row.get("n_resizes")   else 0)
        resize_s.append(float(row["resize_s"])   if row.get("resize_s")    else 0)
        n_ov_res.append(int(row["n_ov_resize"])  if row.get("n_ov_resize") else 0)
        n_load_res.append(int(row["n_load_resize"]) if row.get("n_load_resize") else 0)

m_vals   = np.array(m_vals)
p1       = np.array(p1)
p2       = np.array(p2)
total    = p1 + p2
ov_pct   = np.array(ov_pct)
ov_kmers = np.array(ov_kmers)
n_resizes  = np.array(n_resizes)
resize_s   = np.array(resize_s)
n_ov_res   = np.array(n_ov_res)
n_load_res = np.array(n_load_res)

# Estimated resize wall-time: resize work is single-threaded per partition,
# but partitions run in parallel across t=4 threads → divide by 4.
N_THREADS = 4
resize_wall = resize_s / N_THREADS
non_resize_p2 = np.maximum(p2 - resize_wall, 0)

# ── Figure: 2×2 grid ──────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(title, fontsize=13, fontweight="bold")

x = np.arange(len(m_vals))
xlabels = [f"m={v}" for v in m_vals]
w = 0.55

# colour palette
C_P1     = "#4e79a7"   # blue   — phase1
C_P2     = "#f28e2b"   # orange — phase2
C_OV     = "#e15759"   # red    — overflow
C_NORES  = "#76b7b2"   # teal   — non-resize phase2
C_RWALL  = "#ff9da7"   # pink   — resize wall time
C_OVR    = "#e15759"   # red    — overflow-triggered resizes
C_LOADR  = "#b07aa1"   # purple — load-triggered resizes
C_RS     = "#59a14f"   # green  — resize CPU time line

# ── TL: Phase timing breakdown ────────────────────────────────────────────────
ax = axes[0, 0]
ax.bar(x, p1, width=w, label="phase 1 (partition)", color=C_P1)
ax.bar(x, p2, width=w, bottom=p1, label="phase 2 (count+write)", color=C_P2)
for xi, tot in zip(x, total):
    ax.text(xi, tot + 3, f"{tot:.0f}s", ha="center", va="bottom", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(xlabels)
ax.set_ylabel("Wall time (s)")
ax.set_title("Phase timing breakdown")
ax.legend(loc="upper right", fontsize=9)
if kmc_wall is not None:
    ax.axhline(kmc_wall, color="#e15759", linestyle="--", linewidth=1.8, label=f"KMC {kmc_wall:.0f}s")
    ax.text(len(x) - 0.5, kmc_wall + max(total) * 0.015, f"KMC  {kmc_wall:.0f}s",
            ha="right", va="bottom", fontsize=9, color="#e15759")
    ax.legend(loc="upper left", fontsize=9)
ax.set_ylim(0, max(max(total), kmc_wall if kmc_wall else 0) * 1.18)
ax.yaxis.set_minor_locator(ticker.MultipleLocator(20))
ax.grid(axis="y", alpha=0.3)

# ── TR: Resize events (stacked ov/load) + resize CPU cost ────────────────────
ax = axes[0, 1]
ax.bar(x, n_ov_res,   width=w, label="overflow-triggered", color=C_OVR, alpha=0.9)
ax.bar(x, n_load_res, width=w, bottom=n_ov_res, label="load-triggered", color=C_LOADR, alpha=0.9)
for xi, nr in zip(x, n_resizes):
    ax.text(xi, nr + max(n_resizes) * 0.01, f"{nr}", ha="center", va="bottom", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(xlabels)
ax.set_ylabel("Resize events (count)", color="black")
ax.set_title("Resize events & CPU cost")
ax.legend(loc="upper right", fontsize=9)
ax.set_ylim(0, max(n_resizes) * 1.18)
ax.grid(axis="y", alpha=0.3)

axr = ax.twinx()
axr.plot(x, resize_s, "s--", color=C_RS, linewidth=2, markersize=7, label="resize CPU (s)")
for xi, rs in zip(x, resize_s):
    axr.text(xi, rs + max(resize_s) * 0.02, f"{rs:.0f}s",
             ha="center", va="bottom", fontsize=8, color=C_RS)
axr.set_ylabel("Total resize CPU time (s)", color=C_RS)
axr.tick_params(axis="y", labelcolor=C_RS)
axr.set_ylim(0, max(resize_s) * 1.3)
axr.legend(loc="upper center", fontsize=9)

# ── BL: Overflow rate + phase2 correlation ───────────────────────────────────
ax = axes[1, 0]
bars = ax.bar(x, ov_pct, width=w, color=C_OV, alpha=0.85, label="overflow %")
for xi, pct, cnt in zip(x, ov_pct, ov_kmers):
    label = f"{cnt/1e6:.0f}M\n({pct:.1f}%)"
    ax.text(xi, pct + max(ov_pct) * 0.02, label,
            ha="center", va="bottom", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(xlabels)
ax.set_ylabel("Overflow k-mers (%)", color=C_OV)
ax.tick_params(axis="y", labelcolor=C_OV)
ax.set_title("Overflow rate vs phase 2 time")
ax.set_ylim(0, max(ov_pct) * 1.3)
ax.grid(axis="y", alpha=0.3)

axr = ax.twinx()
axr.plot(x, p2, "o--", color=C_P2, linewidth=2, markersize=6, label="phase2 (s)")
axr.set_ylabel("phase 2 time (s)", color=C_P2)
axr.tick_params(axis="y", labelcolor=C_P2)
axr.set_ylim(0, max(p2) * 1.3)
axr.legend(loc="upper right", fontsize=9)

# ── BR: Phase2 decomposition: resize wall-est vs non-resize ──────────────────
ax = axes[1, 1]
ax.bar(x, non_resize_p2, width=w, label="phase2 non-resize (est.)", color=C_NORES)
ax.bar(x, resize_wall,   width=w, bottom=non_resize_p2,
       label=f"resize wall (resize_cpu / {N_THREADS})", color=C_RWALL)
for xi, rw, tot in zip(x, resize_wall, p2):
    pct = 100 * rw / tot if tot > 0 else 0
    ax.text(xi, tot + 2, f"{pct:.0f}%\nof p2", ha="center", va="bottom", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(xlabels)
ax.set_ylabel("Estimated wall time (s)")
ax.set_title(f"Phase 2 decomposition (resize est. = resize_cpu / {N_THREADS})")
ax.legend(loc="upper right", fontsize=9)
ax.set_ylim(0, max(p2) * 1.25)
ax.grid(axis="y", alpha=0.3)

# ── Save ──────────────────────────────────────────────────────────────────────

out = Path(csv_path).parent / "m_sweep.png"
plt.tight_layout()
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
plt.show()
