#!/usr/bin/env python3
"""
plot_human_scale.py — Analyse and plot the human genome scale experiment results.

Data source:
  benchmark/results/human_scale/results_*/bench_scale.csv   — tuna + KMC scaling with N files
  benchmark/results/human_scale/results_*/table_stats_N*.csv — per-partition stats (tuna only)
"""

import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
SCALE_DIR   = sorted(glob.glob(os.path.join(RESULTS_DIR, "human_scale", "results_*")))[-1]
OUT         = os.path.join(SCALE_DIR, "analysis.png")

# ── Load bench_scale.csv ──────────────────────────────────────────────────────

df_all = pd.read_csv(os.path.join(SCALE_DIR, "bench_scale.csv"))
df_all = df_all.sort_values(["tool", "n_files"]).reset_index(drop=True)

has_kmc = "kmc" in df_all["tool"].values

df  = df_all[df_all["tool"] == "tuna"].copy().reset_index(drop=True)
kmc = df_all[df_all["tool"] == "kmc" ].copy().reset_index(drop=True) if has_kmc else None

for c in ["n_parts","dbg_load_mean","dbg_load_min","dbg_load_max",
          "dbg_oversize_mean","dbg_unique_mean","dbg_n_resizes","dbg_parts_resized"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

df["phase2_per_unique_ns"] = df["phase2_s"] / df["unique_kmers"] * 1e9
df["phase2_per_part_ms"]   = df["phase2_s"] / df["n_parts"] * 1e3

# ── Load table_stats CSVs ─────────────────────────────────────────────────────

total_inserted = {}
for csv_path in glob.glob(os.path.join(SCALE_DIR, "table_stats_N*.csv")):
    n = int(os.path.basename(csv_path).replace("table_stats_N","").replace(".csv",""))
    ts = pd.read_csv(csv_path)
    total_inserted[n] = ts["n_inserted"].sum()

df["total_inserted"] = df["n_files"].map(total_inserted)
df["coverage_ratio"] = df["total_inserted"] / df["unique_kmers"]

# ── Figure layout ─────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(18, 12))
fig.suptitle("Tuna vs KMC — Human genome scale experiment  (k=31, m=21, 8 threads)",
             fontsize=14, fontweight="bold", y=0.98)

gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

ax_cmp    = fig.add_subplot(gs[0, 0:2])  # tuna vs KMC wall time — wide
ax_nparts = fig.add_subplot(gs[0, 2])    # n_parts vs n_files
ax_wall   = fig.add_subplot(gs[1, 0])    # tuna wall breakdown
ax_lf     = fig.add_subplot(gs[1, 1])    # load factor
ax_cov    = fig.add_subplot(gs[1, 2])    # coverage ratio
ax_p2eff  = fig.add_subplot(gs[2, 0])    # phase2 per unique k-mer
ax_p2pp   = fig.add_subplot(gs[2, 1])    # phase2 per partition
ax_rss    = fig.add_subplot(gs[2, 2])    # RSS comparison

x = df["n_files"].values

# ── Plot 1 (wide): tuna vs KMC wall time comparison ──────────────────────────
ax_cmp.plot(x, df["wall_s"], "o-", color="#4C9BE8", linewidth=2.5, markersize=7,
            label="tuna (wall)")
ax_cmp.plot(x, df["phase1_s"], "s--", color="#4C9BE8", linewidth=1.2, markersize=5,
            alpha=0.6, label="tuna phase1")
ax_cmp.plot(x, df["phase2_s"], "^--", color="#4C9BE8", linewidth=1.2, markersize=5,
            alpha=0.6, label="tuna phase2")

if kmc is not None and not kmc.empty:
    xk = kmc["n_files"].values
    ax_cmp.plot(xk, kmc["wall_s"], "o-", color="#E8834C", linewidth=2.5, markersize=7,
                label="KMC (wall = count + dump)")
    ax_cmp.plot(xk, kmc["phase1_s"], "s--", color="#E8834C", linewidth=1.2, markersize=5,
                alpha=0.6, label="KMC count")
    ax_cmp.plot(xk, kmc["phase2_s"], "^--", color="#E8834C", linewidth=1.2, markersize=5,
                alpha=0.6, label="KMC dump")

ax_cmp.set_xlabel("N files"); ax_cmp.set_ylabel("time (s)")
ax_cmp.set_title("Tuna vs KMC — wall time scaling")
ax_cmp.legend(fontsize=8, ncol=2); ax_cmp.grid(True, alpha=0.3)
ax_cmp.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ── Plot 2: n_parts vs n_files ────────────────────────────────────────────────
ax_nparts.plot(x, df["n_parts"], "o-", color="#6A3D9A", linewidth=2, markersize=6,
               label="auto-tuned n_parts")
ax_nparts.set_xlabel("N files"); ax_nparts.set_ylabel("n_parts")
ax_nparts.set_title("Auto-tuned partition count")
ax_nparts.legend(fontsize=8); ax_nparts.grid(True, alpha=0.3)
ax_nparts.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax_nparts.yaxis.set_major_formatter(mticker.ScalarFormatter())

# ── Plot 3: tuna wall time phase breakdown ────────────────────────────────────
ax_wall.stackplot(x, df["phase1_s"], df["phase2_s"],
                  labels=["phase1 (partition)", "phase2 (count)"],
                  colors=["#4C9BE8", "#E8834C"], alpha=0.85)
ax_wall.plot(x, df["wall_s"], "k-", linewidth=1.2, label="wall total")
ax_wall.set_xlabel("N files"); ax_wall.set_ylabel("time (s)")
ax_wall.set_title("Tuna wall time breakdown")
ax_wall.legend(fontsize=8); ax_wall.grid(True, alpha=0.3)
ax_wall.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ── Plot 4: Load factor ────────────────────────────────────────────────────────
if df["dbg_load_mean"].notna().any():
    ax_lf.fill_between(x, df["dbg_load_min"], df["dbg_load_max"],
                        alpha=0.25, color="#FF7F00", label="min–max range")
    ax_lf.plot(x, df["dbg_load_mean"], "o-", color="#FF7F00", linewidth=2, markersize=6,
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
ax_lf.set_xlabel("N files"); ax_lf.set_ylabel("load factor")
ax_lf.set_title("Hash table load factor (tuna)")
ax_lf.legend(fontsize=8); ax_lf.grid(True, alpha=0.3)
ax_lf.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ── Plot 5: Coverage ratio ─────────────────────────────────────────────────────
if df["coverage_ratio"].notna().any():
    ax_cov.plot(x, df["coverage_ratio"], "P-", color="#E31A1C", linewidth=2, markersize=7)
    ax_cov.set_xlabel("N files"); ax_cov.set_ylabel("total_inserted / unique_kmers")
    ax_cov.set_title("Avg k-mer multiplicity (coverage ratio)")
    ax_cov.grid(True, alpha=0.3)
    ax_cov.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax_cov.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax_cov.text(0.05, 0.92,
                "Human: unique k-mers saturate fast\n(high within-species conservation)",
                transform=ax_cov.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))
else:
    ax_cov.text(0.5, 0.5, "total_inserted not available\n(table_stats CSVs missing)",
                ha="center", va="center", transform=ax_cov.transAxes, fontsize=10, color="gray")
    ax_cov.set_title("Coverage ratio (data missing)")

# ── Plot 6: Phase2 cost per unique k-mer ─────────────────────────────────────
ax_p2eff.plot(x, df["phase2_per_unique_ns"], "D-", color="#1F78B4", linewidth=2, markersize=6)
ax_p2eff.set_xlabel("N files"); ax_p2eff.set_ylabel("phase2 / unique_kmers (ns)")
ax_p2eff.set_title("Tuna phase 2 cost per unique k-mer")
ax_p2eff.grid(True, alpha=0.3)
ax_p2eff.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax_p2eff.yaxis.set_major_formatter(mticker.ScalarFormatter())
first, last = df.iloc[0], df.iloc[-1]
ratio = last["phase2_per_unique_ns"] / first["phase2_per_unique_ns"]
ax_p2eff.annotate(f"×{ratio:.1f} per unique k-mer",
                   xy=(last["n_files"], last["phase2_per_unique_ns"]),
                   xytext=(-40, -30), textcoords="offset points",
                   fontsize=8, color="darkred",
                   arrowprops=dict(arrowstyle="->", color="darkred"))

# ── Plot 7: Phase2 time per partition ─────────────────────────────────────────
ax_p2pp.plot(x, df["phase2_per_part_ms"], "v-", color="#B15928", linewidth=2, markersize=6)
ax_p2pp.set_xlabel("N files"); ax_p2pp.set_ylabel("phase2_wall / n_parts (ms/partition)")
ax_p2pp.set_title("Tuna phase 2 wall time per partition")
ax_p2pp.grid(True, alpha=0.3)
ax_p2pp.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ── Plot 8: RSS comparison ─────────────────────────────────────────────────────
ax_rss.plot(x, df["rss_mb"] / 1024, "o-", color="#4C9BE8", linewidth=2, markersize=6,
            label="tuna")
if kmc is not None and not kmc.empty:
    xk = kmc["n_files"].values
    ax_rss.plot(xk, kmc["rss_mb"] / 1024, "o-", color="#E8834C", linewidth=2, markersize=6,
                label="KMC")
ax_rss.set_xlabel("N files"); ax_rss.set_ylabel("peak RSS (GB)")
ax_rss.set_title("Peak memory usage")
ax_rss.legend(fontsize=8); ax_rss.grid(True, alpha=0.3)
ax_rss.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ── Save ──────────────────────────────────────────────────────────────────────

plt.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT}")

# ── Print summary table ────────────────────────────────────────────────────────

print("\n=== tuna summary ===")
tuna_cols = ["n_files","n_parts","phase1_s","phase2_s","wall_s","rss_mb","unique_kmers",
             "dbg_load_mean","dbg_oversize_mean","dbg_n_resizes",
             "phase2_per_unique_ns","coverage_ratio"]
print(df[[c for c in tuna_cols if c in df.columns]].to_string(index=False))

if kmc is not None and not kmc.empty:
    print("\n=== KMC summary ===")
    kmc_cols = ["n_files","phase1_s","phase2_s","wall_s","rss_mb","unique_kmers"]
    print(kmc[[c for c in kmc_cols if c in kmc.columns]].to_string(index=False))

    merged = pd.merge(
        df[["n_files","wall_s","rss_mb","phase1_s","phase2_s"]].rename(
            columns={"wall_s":"tuna_wall","rss_mb":"tuna_rss_mb",
                     "phase1_s":"tuna_p1","phase2_s":"tuna_p2"}),
        kmc[["n_files","wall_s","rss_mb","phase1_s","phase2_s"]].rename(
            columns={"wall_s":"kmc_wall","rss_mb":"kmc_rss_mb",
                     "phase1_s":"kmc_p1","phase2_s":"kmc_p2"}),
        on="n_files", how="inner")
    merged["tuna/kmc_wall"] = merged["tuna_wall"] / merged["kmc_wall"]
    merged["tuna/kmc_rss"]  = merged["tuna_rss_mb"] / merged["kmc_rss_mb"]
    print("\n=== tuna vs KMC — wall time and RSS ratios ===")
    print(merged[["n_files","tuna_wall","kmc_wall","tuna/kmc_wall",
                  "tuna_rss_mb","kmc_rss_mb","tuna/kmc_rss"]].to_string(index=False))
