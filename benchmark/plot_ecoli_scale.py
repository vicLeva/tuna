#!/usr/bin/env python3
"""
plot_ecoli_scale.py — Analyse and plot the ecoli scale experiment results.

Data sources:
  benchmark/results/ecoli_scale/results_*/bench_scale.csv   — tuna + KMC scaling with N files
  benchmark/results/ecoli_scale/results_*/table_stats_N*.csv — per-partition stats (tuna only)
  benchmark/results/n_sweep_ecoli.tsv                        — partition count sweep (200 E. coli)
"""

import sys
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
SCALE_DIR   = sorted(glob.glob(os.path.join(RESULTS_DIR, "ecoli_scale", "results_*")))[-1]
OUT         = os.path.join(SCALE_DIR, "analysis.png")

# ── Load bench_scale.csv ──────────────────────────────────────────────────────

df_all = pd.read_csv(os.path.join(SCALE_DIR, "bench_scale.csv"))
df_all = df_all.sort_values(["tool", "n_files"]).reset_index(drop=True)

has_tool_col = "tool" in df_all.columns
has_kmc      = has_tool_col and "kmc" in df_all["tool"].values

if has_tool_col:
    df  = df_all[df_all["tool"] == "tuna"].copy().reset_index(drop=True)
    kmc = df_all[df_all["tool"] == "kmc" ].copy().reset_index(drop=True) if has_kmc else None
else:
    # Legacy: no tool column — all rows are tuna
    df  = df_all.copy()
    kmc = None

# Numeric coerce for na-filled KMC columns
for c in ["n_parts","dbg_load_mean","dbg_load_min","dbg_load_max",
          "dbg_oversize_mean","dbg_unique_mean","dbg_n_resizes","dbg_parts_resized"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# Derived tuna columns
df["phase2_per_unique_ns"] = df["phase2_s"] / df["unique_kmers"] * 1e9
df["phase2_per_part_ms"]   = df["phase2_s"] / df["n_parts"] * 1e3
df["phase1_per_file_s"]    = df["phase1_s"] / df["n_files"]

# ── Load table_stats CSVs → total insertions per N ───────────────────────────

total_inserted = {}
for csv_path in glob.glob(os.path.join(SCALE_DIR, "table_stats_N*.csv")):
    n = int(os.path.basename(csv_path).replace("table_stats_N","").replace(".csv",""))
    ts = pd.read_csv(csv_path)
    total_inserted[n] = ts["n_inserted"].sum()

df["total_inserted"] = df["n_files"].map(total_inserted)
df["coverage_ratio"] = df["total_inserted"] / df["unique_kmers"]

# ── Load n_sweep_ecoli.tsv ────────────────────────────────────────────────────

sweep_path = os.path.join(RESULTS_DIR, "n_sweep_ecoli.tsv")
sw = None
if os.path.exists(sweep_path):
    sw = pd.read_csv(sweep_path, sep="\t")
    sw = sw.groupby("n")[["phase1","phase2","total"]].mean().reset_index()

# ── Figure layout ─────────────────────────────────────────────────────────────
# Row 0: tuna vs KMC wall time comparison (wide) | n_parts vs N
# Row 1: phase2 breakdown | load factor | coverage ratio
# Row 2: phase2/unique k-mer | phase2/partition | n_sweep

fig = plt.figure(figsize=(18, 14))
fig.suptitle("Tuna vs KMC — E. coli scale experiment  (k=31, m=21, 8 threads)",
             fontsize=14, fontweight="bold", y=0.98)

gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

ax_cmp    = fig.add_subplot(gs[0, 0:2])  # tuna vs KMC wall time — wide
ax_nparts = fig.add_subplot(gs[0, 2])    # n_parts vs n_files
ax_wall   = fig.add_subplot(gs[1, 0])    # tuna wall time breakdown (phase1+phase2)
ax_lf     = fig.add_subplot(gs[1, 1])    # load factor
ax_cov    = fig.add_subplot(gs[1, 2])    # coverage ratio
ax_p2eff  = fig.add_subplot(gs[2, 0])    # phase2 per unique k-mer
ax_p2pp   = fig.add_subplot(gs[2, 1])    # phase2 per partition
ax_sweep  = fig.add_subplot(gs[2, 2])    # n_sweep

x = df["n_files"].values

# ── Plot 1 (wide): tuna vs KMC wall time comparison ──────────────────────────
ax_cmp.plot(x, df["wall_s"], "o-", color="#4C9BE8", linewidth=2.5, markersize=6,
            label="tuna (wall)")
ax_cmp.plot(x, df["phase1_s"], "s--", color="#4C9BE8", linewidth=1.2, markersize=4,
            alpha=0.6, label="tuna phase1")
ax_cmp.plot(x, df["phase2_s"], "^--", color="#4C9BE8", linewidth=1.2, markersize=4,
            alpha=0.6, label="tuna phase2")

if kmc is not None and not kmc.empty:
    xk = kmc["n_files"].values
    ax_cmp.plot(xk, kmc["wall_s"], "o-", color="#E8834C", linewidth=2.5, markersize=6,
                label="KMC (wall = count + dump)")
    ax_cmp.plot(xk, kmc["phase1_s"], "s--", color="#E8834C", linewidth=1.2, markersize=4,
                alpha=0.6, label="KMC count")
    ax_cmp.plot(xk, kmc["phase2_s"], "^--", color="#E8834C", linewidth=1.2, markersize=4,
                alpha=0.6, label="KMC dump")

ax_cmp.set_xscale("log"); ax_cmp.set_yscale("log")
ax_cmp.set_xlabel("N files"); ax_cmp.set_ylabel("time (s)")
ax_cmp.set_title("Tuna vs KMC — wall time scaling")
ax_cmp.legend(fontsize=8, ncol=2); ax_cmp.grid(True, alpha=0.3)
ax_cmp.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 2: n_parts vs n_files ────────────────────────────────────────────────
ax_nparts.plot(x, df["n_parts"], "o-", color="#6A3D9A", linewidth=2, markersize=5,
               label="auto-tuned n_parts")
y0 = df["n_parts"].iloc[0]; x0 = x[0]
ax_nparts.plot(x, y0 * x / x0, "k--", linewidth=1, alpha=0.5, label="O(N) reference")
ax_nparts.set_xscale("log"); ax_nparts.set_yscale("log")
ax_nparts.set_xlabel("N files"); ax_nparts.set_ylabel("n_parts")
ax_nparts.set_title("Auto-tuned partition count")
ax_nparts.legend(fontsize=8); ax_nparts.grid(True, alpha=0.3)
ax_nparts.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax_nparts.yaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 3: tuna wall time phase breakdown ────────────────────────────────────
ax_wall.stackplot(x, df["phase1_s"], df["phase2_s"],
                  labels=["phase1 (partition)", "phase2 (count)"],
                  colors=["#4C9BE8", "#E8834C"], alpha=0.85)
ax_wall.plot(x, df["wall_s"], "k-", linewidth=1.2, label="wall total")
ax_wall.set_xscale("log"); ax_wall.set_yscale("log")
ax_wall.set_xlabel("N files"); ax_wall.set_ylabel("time (s)")
ax_wall.set_title("Tuna wall time breakdown")
ax_wall.legend(fontsize=8); ax_wall.grid(True, alpha=0.3)
ax_wall.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 4: Load factor ────────────────────────────────────────────────────────
if df["dbg_load_mean"].notna().any():
    ax_lf.fill_between(x, df["dbg_load_min"], df["dbg_load_max"],
                        alpha=0.25, color="#FF7F00", label="min–max range")
    ax_lf.plot(x, df["dbg_load_mean"], "o-", color="#FF7F00", linewidth=2, markersize=5,
               label="mean load factor")
    ax_lf.axhline(0.8, color="red", linestyle="--", linewidth=1, alpha=0.6, label="0.8 resize trigger")
    ax_lf.set_ylim(0, 1.05)
    for _, row in df[df["dbg_n_resizes"] > 0].iterrows():
        ax_lf.annotate("resize!", xy=(row["n_files"], row["dbg_load_max"]),
                        xytext=(0, 8), textcoords="offset points",
                        ha="center", fontsize=7, color="red")
else:
    ax_lf.text(0.5, 0.5, "no dbg data", ha="center", va="center",
               transform=ax_lf.transAxes, color="gray")
ax_lf.set_xscale("log")
ax_lf.set_xlabel("N files"); ax_lf.set_ylabel("load factor")
ax_lf.set_title("Hash table load factor (tuna)")
ax_lf.legend(fontsize=8); ax_lf.grid(True, alpha=0.3)
ax_lf.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 5: Coverage ratio ─────────────────────────────────────────────────────
if df["coverage_ratio"].notna().any():
    ax_cov.plot(x, df["coverage_ratio"], "P-", color="#E31A1C", linewidth=2, markersize=6)
    ax_cov.set_xscale("log"); ax_cov.set_yscale("log")
    ax_cov.set_xlabel("N files"); ax_cov.set_ylabel("total_inserted / unique_kmers")
    ax_cov.set_title("Avg k-mer multiplicity (coverage ratio)")
    ax_cov.grid(True, alpha=0.3)
    ax_cov.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax_cov.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax_cov.text(0.05, 0.92,
                "High multiplicity → Phase 2\nre-probes same slots many times",
                transform=ax_cov.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))
else:
    ax_cov.text(0.5, 0.5, "total_inserted not available\n(table_stats CSVs missing)",
                ha="center", va="center", transform=ax_cov.transAxes, fontsize=10, color="gray")
    ax_cov.set_title("Coverage ratio (data missing)")

# ── Plot 6: Phase2 efficiency — cost per unique k-mer ────────────────────────
ax_p2eff.plot(x, df["phase2_per_unique_ns"], "D-", color="#1F78B4", linewidth=2, markersize=5)
ax_p2eff.set_xscale("log"); ax_p2eff.set_yscale("log")
ax_p2eff.set_xlabel("N files"); ax_p2eff.set_ylabel("phase2 / unique_kmers (ns)")
ax_p2eff.set_title("Tuna phase 2 cost per unique k-mer")
ax_p2eff.grid(True, alpha=0.3)
ax_p2eff.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax_p2eff.yaxis.set_major_formatter(mticker.ScalarFormatter())
first, last = df.iloc[0], df.iloc[-1]
ratio = last["phase2_per_unique_ns"] / first["phase2_per_unique_ns"]
ax_p2eff.annotate(f"×{ratio:.0f} per unique k-mer",
                   xy=(last["n_files"], last["phase2_per_unique_ns"]),
                   xytext=(-60, -30), textcoords="offset points",
                   fontsize=8, color="darkred",
                   arrowprops=dict(arrowstyle="->", color="darkred"))

# ── Plot 7: Phase2 time per partition ─────────────────────────────────────────
ax_p2pp.plot(x, df["phase2_per_part_ms"], "v-", color="#B15928", linewidth=2, markersize=5)
ax_p2pp.set_xscale("log")
ax_p2pp.set_xlabel("N files"); ax_p2pp.set_ylabel("phase2_wall / n_parts (ms/partition)")
ax_p2pp.set_title("Tuna phase 2 wall time per partition")
ax_p2pp.grid(True, alpha=0.3)
ax_p2pp.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 8: n_sweep — phase2 vs n_parts ───────────────────────────────────────
if sw is not None:
    ax_sweep.plot(sw["n"], sw["phase2"], "o-", color="#6A3D9A", linewidth=2, markersize=6, label="phase2")
    ax_sweep.plot(sw["n"], sw["phase1"], "s--", color="#4C9BE8", linewidth=1.5, markersize=5, label="phase1")
    ax_sweep.plot(sw["n"], sw["total"], "k-", linewidth=1.5, alpha=0.5, label="total")
    n200_nparts = df[df["n_files"] == 200]["n_parts"].values
    if len(n200_nparts):
        ax_sweep.axvline(n200_nparts[0], color="red", linestyle=":", linewidth=1.5,
                          label=f"auto n={n200_nparts[0]} (N=200 files)")
    ax_sweep.set_xscale("log")
    ax_sweep.set_xlabel("n_parts"); ax_sweep.set_ylabel("time (s)")
    ax_sweep.set_title("Partition count sweep — fixed 200 E. coli files")
    ax_sweep.legend(fontsize=8); ax_sweep.grid(True, alpha=0.3)
    ax_sweep.xaxis.set_major_formatter(mticker.ScalarFormatter())
else:
    ax_sweep.text(0.5, 0.5, "n_sweep_ecoli.tsv not found", ha="center", va="center",
                  transform=ax_sweep.transAxes, fontsize=11, color="gray")
    ax_sweep.set_title("Partition count sweep (data missing)")

# ── Save ──────────────────────────────────────────────────────────────────────

plt.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT}")

# ── Print summary table ────────────────────────────────────────────────────────

print("\n=== tuna summary ===")
tuna_cols = ["n_files","n_parts","phase1_s","phase2_s","wall_s","unique_kmers",
             "dbg_load_mean","dbg_oversize_mean","dbg_n_resizes",
             "phase2_per_unique_ns","coverage_ratio"]
print(df[[c for c in tuna_cols if c in df.columns]].to_string(index=False))

if kmc is not None and not kmc.empty:
    print("\n=== KMC summary ===")
    kmc_cols = ["n_files","phase1_s","phase2_s","wall_s","unique_kmers"]
    print(kmc[[c for c in kmc_cols if c in kmc.columns]].to_string(index=False))

    # Side-by-side comparison on shared N values
    merged = pd.merge(
        df[["n_files","wall_s","phase1_s","phase2_s"]].rename(
            columns={"wall_s":"tuna_wall","phase1_s":"tuna_p1","phase2_s":"tuna_p2"}),
        kmc[["n_files","wall_s","phase1_s","phase2_s"]].rename(
            columns={"wall_s":"kmc_wall","phase1_s":"kmc_p1","phase2_s":"kmc_p2"}),
        on="n_files", how="inner")
    merged["tuna/kmc"] = merged["tuna_wall"] / merged["kmc_wall"]
    print("\n=== tuna vs KMC — wall time ratio ===")
    print(merged[["n_files","tuna_wall","kmc_wall","tuna/kmc"]].to_string(index=False))
