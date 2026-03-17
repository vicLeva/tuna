#!/usr/bin/env python3
"""Plot E. coli N-sweep benchmark: tuna vs KMC wall time + insert/unique ratio."""

import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# ── Load data ─────────────────────────────────────────────────────────────────

csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    sorted(Path.home().glob("tuna_bench_tmp/bench_n_sweep_*/bench_n.csv"))[-1]
title = sys.argv[2] if len(sys.argv) > 2 else \
    "tuna vs KMC — E. coli genomes, k=31, m=17, t=4"

n_files = []
tuna_p1, tuna_p2, tuna_total = [], [], []
tuna_total_k, tuna_unique_k = [], []
kmc_count_s, kmc_dump_s, kmc_total = [], [], []
kmc_unique_k = []

with open(csv_path) as f:
    for row in csv.DictReader(f):
        n_files.append(int(row["n_files"]))
        tuna_p1.append(float(row["tuna_p1_s"]))
        tuna_p2.append(float(row["tuna_p2_s"]))
        tuna_total.append(float(row["tuna_total_s"]))
        tuna_total_k.append(int(row["tuna_total_kmers"]) if row["tuna_total_kmers"] else 0)
        tuna_unique_k.append(int(row["tuna_unique_kmers"]) if row["tuna_unique_kmers"] else 0)
        kmc_count_s.append(float(row["kmc_count_s"]))
        kmc_dump_s.append(float(row["kmc_dump_s"]))
        kmc_total.append(float(row["kmc_total_s"]))
        kmc_unique_k.append(int(row["kmc_unique_kmers"]) if row["kmc_unique_kmers"] else 0)

n_files     = np.array(n_files,     dtype=float)
tuna_p1     = np.array(tuna_p1)
tuna_p2     = np.array(tuna_p2)
tuna_total  = np.array(tuna_total)
tuna_total_k  = np.array(tuna_total_k,  dtype=float)
tuna_unique_k = np.array(tuna_unique_k, dtype=float)
kmc_count_s = np.array(kmc_count_s)
kmc_dump_s  = np.array(kmc_dump_s)
kmc_total   = np.array(kmc_total)
kmc_unique_k  = np.array(kmc_unique_k,  dtype=float)

ratio = np.where(tuna_unique_k > 0, tuna_total_k / tuna_unique_k, np.nan)
speedup = kmc_total / tuna_total   # >1 → tuna faster

# ── Colour palette ────────────────────────────────────────────────────────────
C_P1      = "#4e79a7"   # blue   — tuna phase1
C_P2      = "#76b7b2"   # teal   — tuna phase2
C_TUNA    = "#f28e2b"   # orange — tuna total
C_KMC_CNT = "#b07aa1"   # purple — KMC count
C_KMC_DMP = "#ff9da7"   # pink   — KMC dump
C_KMC     = "#e15759"   # red    — KMC total
C_RATIO   = "#59a14f"   # green  — ratio / unique

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(title, fontsize=13, fontweight="bold")

xlabels = [str(int(n)) for n in n_files]

# ── TL: Wall time — log-log, tuna p1/p2/total + KMC count/dump/total ─────────
ax = axes[0, 0]

# Phase breakdown (dashed, thinner)
ax.plot(n_files, tuna_p1,     "o--", color=C_P1,      linewidth=1.5, markersize=5,
        alpha=0.8, label="tuna phase1 (partition)")
ax.plot(n_files, tuna_p2,     "s--", color=C_P2,      linewidth=1.5, markersize=5,
        alpha=0.8, label="tuna phase2 (count+write)")
ax.plot(n_files, kmc_count_s, "^--", color=C_KMC_CNT, linewidth=1.5, markersize=5,
        alpha=0.8, label="KMC count")
ax.plot(n_files, kmc_dump_s,  "v--", color=C_KMC_DMP, linewidth=1.5, markersize=5,
        alpha=0.8, label="KMC dump")

# Totals (solid, thick)
ax.plot(n_files, tuna_total,  "o-",  color=C_TUNA, linewidth=2.5, markersize=8,
        label="tuna total", zorder=5)
ax.plot(n_files, kmc_total,   "s-",  color=C_KMC,  linewidth=2.5, markersize=8,
        label="KMC total",  zorder=5)

# Labels on totals at each point
for n, tt, kt in zip(n_files, tuna_total, kmc_total):
    ax.annotate(f"{tt:.0f}s", (n, tt), textcoords="offset points",
                xytext=(0, 6), ha="center", fontsize=7, color=C_TUNA)
    ax.annotate(f"{kt:.0f}s", (n, kt), textcoords="offset points",
                xytext=(0, -12), ha="center", fontsize=7, color=C_KMC)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(n_files)
ax.set_xticklabels(xlabels)
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xlabel("N genomes (log scale)")
ax.set_ylabel("Wall time, s (log scale)")
ax.set_title("Wall time scaling — log/log")
ax.legend(loc="upper left", fontsize=7, ncol=2)
ax.grid(which="both", alpha=0.25)

# ── TR: Speedup ratio KMC÷tuna — linear y, log x ─────────────────────────────
ax = axes[0, 1]
ax.axhline(1.0, color="black", linestyle="--", linewidth=1.5)
ax.plot(n_files, speedup, "D-", color="#59a14f", linewidth=2, markersize=8)
for n, s in zip(n_files, speedup):
    va = "bottom" if s >= 1 else "top"
    ofs = 4 if s >= 1 else -4
    ax.annotate(f"{s:.2f}×", (n, s), textcoords="offset points",
                xytext=(0, ofs), ha="center", va=va, fontsize=8,
                color=C_TUNA if s > 1 else C_KMC, fontweight="bold")

# Shade regions
ax.fill_between(n_files, 1, speedup,
                where=speedup >= 1, alpha=0.12, color=C_TUNA, label="tuna faster")
ax.fill_between(n_files, speedup, 1,
                where=speedup < 1,  alpha=0.12, color=C_KMC,  label="KMC faster")

ax.set_xscale("log")
ax.set_xticks(n_files)
ax.set_xticklabels(xlabels)
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xlabel("N genomes (log scale)")
ax.set_ylabel("KMC total / tuna total")
ax.set_title("Speedup (KMC÷tuna)  >1 = tuna wins")
ax.legend(loc="upper right", fontsize=8)
ax.grid(which="both", alpha=0.25)

# ── BL: k-mer counts — log-log ────────────────────────────────────────────────
ax = axes[1, 0]
ax.plot(n_files, tuna_total_k  / 1e6, "s--", color=C_P1,   linewidth=2, markersize=7,
        label="total insertions (tuna)")
ax.plot(n_files, tuna_unique_k / 1e6, "o-",  color=C_TUNA, linewidth=2, markersize=7,
        label="unique k-mers (tuna)")
ax.plot(n_files, kmc_unique_k  / 1e6, "^:",  color=C_KMC,  linewidth=1.5, markersize=6,
        label="unique k-mers (KMC)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(n_files)
ax.set_xticklabels(xlabels)
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xlabel("N genomes (log scale)")
ax.set_ylabel("k-mers, millions (log scale)")
ax.set_title("k-mer counts — log/log")
ax.legend(loc="upper left", fontsize=8)
ax.grid(which="both", alpha=0.25)

# ── BR: Redundancy — log x, linear y ─────────────────────────────────────────
ax = axes[1, 1]
ax.plot(n_files, ratio, "o-", color=C_RATIO, linewidth=2, markersize=8)
for n, r in zip(n_files, ratio):
    if not np.isnan(r):
        ax.annotate(f"{r:.0f}×", (n, r), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=8)
ax.set_xscale("log")
ax.set_xticks(n_files)
ax.set_xticklabels(xlabels)
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xlabel("N genomes (log scale)")
ax.set_ylabel("total insertions / unique k-mers")
ax.set_title("k-mer redundancy (insert/unique ratio)")
ax.grid(which="both", alpha=0.25)

# ── Save ──────────────────────────────────────────────────────────────────────
out = Path(csv_path).parent / "ecoli_n_sweep.png"
plt.tight_layout()
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
plt.show()
