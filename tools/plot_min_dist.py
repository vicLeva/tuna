#!/usr/bin/env python3
"""
plot_min_dist.py — plot the minimizer coverage distribution from tuna -dbg output.

Usage:
    python3 tools/plot_min_dist.py <debug_min_coverage.csv> [output.png]

The CSV is written by tuna -dbg to <work_dir>/debug_min_coverage.csv.
Columns: coverage, n_minimizers, total_kmers
  coverage     — number of k-mers sharing one minimizer
  n_minimizers — how many minimizers have that coverage
  total_kmers  — coverage * n_minimizers (sanity check)
"""

import sys
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

def load_csv(path):
    coverage, n_min, total_km = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            coverage.append(int(row["coverage"]))
            n_min.append(int(row["n_minimizers"]))
            total_km.append(int(row["total_kmers"]))
    return np.array(coverage), np.array(n_min), np.array(total_km)

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    csv_path  = sys.argv[1]
    out_path  = sys.argv[2] if len(sys.argv) > 2 else csv_path.replace(".csv", ".png")

    cov, n_min, total_km = load_csv(csv_path)

    total_minimizers = n_min.sum()
    total_kmers      = total_km.sum()
    avg_cov          = total_kmers / total_minimizers if total_minimizers else 0
    max_cov          = cov.max()

    # Cumulative k-mers covered as a function of minimizer coverage threshold.
    # i.e. "what fraction of k-mers come from minimizers with coverage ≤ X?"
    sort_idx   = np.argsort(cov)
    cov_s      = cov[sort_idx]
    km_s       = total_km[sort_idx]
    cum_kmers  = np.cumsum(km_s) / total_kmers * 100.0

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    title_base = os.path.basename(csv_path)
    fig.suptitle(
        f"{title_base}  |  {total_minimizers:,} minimizers  |  {total_kmers:,} k-mers  |  avg={avg_cov:.1f}  max={max_cov:,}",
        fontsize=10
    )

    # ── Left: histogram of coverage (log-log) ──────────────────────────────────
    ax = axes[0]
    # Bin by powers of 2 for a log-histogram.
    bins = np.unique(np.concatenate([[1], 2**np.arange(0, int(np.log2(max_cov)) + 2)]))
    bin_n_min  = np.zeros(len(bins) - 1, dtype=np.int64)
    bin_n_km   = np.zeros(len(bins) - 1, dtype=np.int64)
    for i in range(len(bins) - 1):
        mask = (cov >= bins[i]) & (cov < bins[i + 1])
        bin_n_min[i] = n_min[mask].sum()
        bin_n_km[i]  = total_km[mask].sum()

    bin_centers = np.sqrt(bins[:-1] * bins[1:])
    bin_widths  = bins[1:] - bins[:-1]
    ax.bar(bin_centers, bin_n_min, width=bin_widths * 0.8, align="center",
           color="steelblue", alpha=0.8, label="# minimizers")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k-mers per minimizer (coverage)", fontsize=11)
    ax.set_ylabel("number of minimizers", fontsize=11)
    ax.set_title("Minimizer coverage distribution (log-log)", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation(base=10, labelOnlyBase=False))
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend()

    # ── Right: cumulative k-mers vs coverage threshold ─────────────────────────
    ax2 = axes[1]
    ax2.plot(cov_s, cum_kmers, color="darkorange", linewidth=1.5)
    ax2.set_xscale("log")
    ax2.set_xlabel("minimizer coverage (k-mers per minimizer)", fontsize=11)
    ax2.set_ylabel("cumulative % of all k-mers", fontsize=11)
    ax2.set_title("Cumulative k-mers by minimizer coverage", fontsize=11)
    ax2.set_ylim(0, 100)
    ax2.axhline(50,  color="gray", linestyle="--", linewidth=0.8, label="50%")
    ax2.axhline(90,  color="gray", linestyle=":",  linewidth=0.8, label="90%")
    ax2.axhline(99,  color="red",  linestyle=":",  linewidth=0.8, label="99%")
    ax2.grid(True, which="both", linestyle="--", alpha=0.4)
    ax2.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")

    # ── Text summary ────────────────────────────────────────────────────────────
    # Find coverage values at 50/90/99% of k-mers.
    for pct in (50, 90, 99):
        idx = np.searchsorted(cum_kmers, pct)
        if idx < len(cov_s):
            print(f"  {pct}% of k-mers come from minimizers with coverage ≤ {cov_s[idx]:,}")

if __name__ == "__main__":
    main()
