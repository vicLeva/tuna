#!/usr/bin/env python3
"""
plot_part_modes.py  —  visualise tuna partition-strategy benchmark results.

Reads two CSV files produced by bench_part_modes.sh:
  <out_dir>/runs.csv        — one row per (mode, threads, rep) run
  <out_dir>/partitions.csv  — one row per (run_id, partition_id), rep=1 only

Produces five figures saved under <out_dir>/:
  plot_wall_time.png   — wall time vs thread count, one line per mode
  plot_speedup.png     — parallel speedup relative to 1-thread baseline
  plot_phases.png      — stacked phase-time bars (scan+partition / count+write)
  plot_balance.png     — partition load-balance (max/mean) and per-partition
                         byte-size violin plots
  plot_summary.png     — 2×2 summary grid (wall time, speedup, phases, balance)

Usage:
  python3 plot_part_modes.py <out_dir>
  python3 plot_part_modes.py <out_dir> --dpi 150 --format pdf

Requirements:
  pip install pandas matplotlib numpy
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

# ── Aesthetics ────────────────────────────────────────────────────────────────

PALETTE = {
    "kmc":      "#388E3C",   # green
    "kmtricks": "#F57C00",   # orange
    "hash":     "#1976D2",   # blue
    "xxhash":   "#7B1FA2",   # purple
}
MARKERS = {"kmc": "^", "kmtricks": "s", "hash": "o", "xxhash": "D"}
LABELS  = {
    "kmc":      "KMC (default)",
    "kmtricks": "kmtricks",
    "hash":     "ntHash",
    "xxhash":   "XXH3",
}
MODE_ORDER = ["kmc", "kmtricks", "hash", "xxhash"]

PHASE_COLORS = {
    "scan+partition": "#90CAF9",   # light blue
    "count+write":   "#A5D6A7",   # light green
}


def style_ax(ax, *, xlabel="", ylabel="", title="", grid_axis="both"):
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.grid(True, alpha=0.3, axis=grid_axis, linestyle="--", linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── Data loading & aggregation ────────────────────────────────────────────────

def load_data(out_dir: Path):
    runs_csv  = out_dir / "runs.csv"
    parts_csv = out_dir / "partitions.csv"

    for f in (runs_csv, parts_csv):
        if not f.exists():
            sys.exit(f"ERROR: file not found: {f}")

    runs  = pd.read_csv(runs_csv)
    parts = pd.read_csv(parts_csv)

    # Only keep modes that are present in the data (future-proof for partial runs)
    present = [m for m in MODE_ORDER if m in runs["mode"].values]
    runs["mode"]  = pd.Categorical(runs["mode"],  categories=present, ordered=True)
    parts["mode"] = pd.Categorical(parts["mode"], categories=present, ordered=True)

    # Aggregate over reps: mean + std of numeric columns
    num_cols = [
        "wall_s", "peak_rss_kb",
        "phase01_s", "phase2_s", "total_s",
        "kmers_unique",
        "p_mean_bytes", "p_std_bytes", "p_max_bytes", "p_min_bytes", "p_imbalance",
    ]
    means = (runs.groupby(["mode", "threads"], observed=True)[num_cols]
                 .mean().reset_index())
    stds  = (runs.groupby(["mode", "threads"], observed=True)[num_cols]
                 .std().fillna(0).reset_index())

    all_threads = sorted(runs["threads"].unique())
    min_t = min(all_threads)
    max_t = max(all_threads)

    return runs, parts, means, stds, all_threads, min_t, max_t, present


# ── Plot 1: Wall time vs threads ──────────────────────────────────────────────

def plot_wall_time(ax, means, stds, all_threads, modes=None):
    for mode in (modes or MODE_ORDER):
        dm = means[means["mode"] == mode].sort_values("threads")
        ds = stds [stds ["mode"] == mode].sort_values("threads")
        ax.errorbar(
            dm["threads"], dm["wall_s"],
            yerr=ds["wall_s"],
            color=PALETTE[mode], marker=MARKERS[mode],
            label=LABELS[mode],
            linewidth=2, markersize=7, capsize=4, capthick=1.5,
        )

    ax.set_xticks(all_threads)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
    style_ax(ax,
             xlabel="Threads",
             ylabel="Wall time (s)",
             title="Total wall time vs thread count")
    ax.legend(fontsize=9)


# ── Plot 2: Speedup ───────────────────────────────────────────────────────────

def plot_speedup(ax, means, all_threads, min_t, modes=None):
    ideal_x = all_threads
    ideal_y = [t / min_t for t in ideal_x]
    ax.plot(ideal_x, ideal_y, "k--", alpha=0.35, linewidth=1.5,
            label="Ideal (linear)")

    for mode in (modes or MODE_ORDER):
        dm = means[means["mode"] == mode].sort_values("threads")
        baseline_rows = dm[dm["threads"] == min_t]["wall_s"].values
        if len(baseline_rows) == 0:
            continue
        baseline = baseline_rows[0]
        speedup  = baseline / dm["wall_s"].values
        ax.plot(
            dm["threads"], speedup,
            color=PALETTE[mode], marker=MARKERS[mode],
            label=LABELS[mode],
            linewidth=2, markersize=7,
        )

    ax.set_xticks(all_threads)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
    style_ax(ax,
             xlabel="Threads",
             ylabel=f"Speedup  (relative to {min_t} thread{'s' if min_t > 1 else ''})",
             title="Parallel speedup")
    ax.legend(fontsize=9)


# ── Plot 3: Phase breakdown (one panel per mode) ──────────────────────────────
#
# axes: sequence of Axes objects, one per mode in MODE_ORDER.
# All panels share the same y-scale so bar heights are directly comparable.

def plot_phases(axes, means, all_threads, modes=None):
    modes = modes or MODE_ORDER
    # Shared y-limit: tallest total bar across all modes
    y_max = 0.0
    for mode in modes:
        dm = means[means["mode"] == mode]
        if not dm.empty:
            y_max = max(y_max, (dm["phase01_s"] + dm["phase2_s"]).max())
    y_max *= 1.15

    x     = np.arange(len(all_threads))
    bar_w = 0.55

    for ax, mode in zip(axes, modes):
        dm  = means[means["mode"] == mode].sort_values("threads")
        p01 = dm["phase01_s"].values
        p2  = dm["phase2_s"].values

        ax.bar(x, p01, bar_w,
               color=PHASE_COLORS["scan+partition"],
               edgecolor=PALETTE[mode], linewidth=1.2, alpha=0.85,
               label="Scan + Partition")
        ax.bar(x, p2, bar_w, bottom=p01,
               color=PHASE_COLORS["count+write"],
               edgecolor=PALETTE[mode], linewidth=1.2, alpha=0.85,
               label="Count + Write")

        ax.set_xticks(x)
        ax.set_xticklabels([f"t={t}" for t in all_threads])
        ax.set_ylim(0, y_max)
        style_ax(ax, xlabel="Threads", title=LABELS[mode], grid_axis="y")

    # y-label and legend only on the leftmost panel
    axes[0].set_ylabel("Time (s)", fontsize=10)
    for ax in list(axes)[1:]:
        ax.set_ylabel("")
        ax.tick_params(labelleft=False)

    handles, labels_leg = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels_leg, fontsize=9, loc="upper right")


# ── Plot 4: Partition load balance ───────────────────────────────────────────

def plot_balance(ax, parts, means, modes=None):
    """
    Left y-axis: violin plots of per-partition byte sizes (from partitions.csv).
    Right y-axis: max/mean imbalance ratio (from runs.csv means, first threads entry).
    """
    modes = modes or MODE_ORDER
    # Use mean imbalance aggregated over all thread counts (balance is nearly
    # thread-independent for Phase 1).
    imb_by_mode = (means.groupby("mode", observed=True)["p_imbalance"]
                        .mean().reindex(modes))

    x = np.arange(len(modes))
    width = 0.55

    # Violin plots of per-partition sizes ─────────────────────────────────────
    for mi, mode in enumerate(modes):
        data = parts[parts["mode"] == mode]["size_bytes"].values / 1024  # kB
        if len(data) == 0:
            continue
        vp = ax.violinplot(
            data,
            positions=[x[mi]],
            widths=width * 0.9,
            showmedians=True,
            showextrema=True,
        )
        for body in vp["bodies"]:
            body.set_facecolor(PALETTE[mode])
            body.set_alpha(0.5)
            body.set_edgecolor(PALETTE[mode])
        for part_name in ("cmedians", "cmins", "cmaxes", "cbars"):
            if part_name in vp:
                vp[part_name].set_color(PALETTE[mode])
                vp[part_name].set_linewidth(1.5)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in modes], fontsize=9)
    style_ax(ax,
             ylabel="Partition size (kB)",
             title="Partition load balance",
             grid_axis="y")

    # Imbalance ratio on secondary y-axis ─────────────────────────────────────
    ax2 = ax.twinx()
    ax2.plot(x, imb_by_mode.values,
             "D--", color="#B71C1C", linewidth=1.5, markersize=8,
             label="max/mean ratio", zorder=5)
    ax2.axhline(1.0, color="#B71C1C", linestyle=":", alpha=0.4, linewidth=1)
    ax2.set_ylabel("Imbalance  (max / mean)", fontsize=10, color="#B71C1C")
    ax2.tick_params(axis="y", colors="#B71C1C")
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="upper right", fontsize=9)


# ── Plot 5: Peak RSS ──────────────────────────────────────────────────────────

def plot_rss(ax, means, all_threads, modes=None):
    for mode in (modes or MODE_ORDER):
        dm = means[means["mode"] == mode].sort_values("threads")
        ax.plot(
            dm["threads"], dm["peak_rss_kb"] / 1024,   # → MB
            color=PALETTE[mode], marker=MARKERS[mode],
            label=LABELS[mode], linewidth=2, markersize=7,
        )
    ax.set_xticks(all_threads)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
    style_ax(ax,
             xlabel="Threads",
             ylabel="Peak RSS (MB)",
             title="Peak resident memory vs thread count")
    ax.legend(fontsize=9)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot tuna partition-strategy benchmark results."
    )
    parser.add_argument("out_dir", help="Directory produced by bench_part_modes.sh")
    parser.add_argument("--dpi",    type=int,  default=120, help="Output DPI [120]")
    parser.add_argument("--format", default="png",
                        help="Output format: png, pdf, svg [png]")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    fmt     = args.format
    dpi     = args.dpi

    runs, parts, means, stds, all_threads, min_t, max_t, modes = load_data(out_dir)

    print(f"Loaded {len(runs)} run(s) across "
          f"{runs['mode'].nunique()} mode(s), "
          f"{runs['threads'].nunique()} thread config(s), "
          f"{runs['rep'].max()} rep(s).")
    print(f"Partition data: {len(parts)} rows")

    # ── Individual plots ──────────────────────────────────────────────────────

    def save(fig, name):
        p = out_dir / f"{name}.{fmt}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        print(f"  Saved: {p}")
        plt.close(fig)

    n_modes = len(modes)

    # 1. Wall time
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_wall_time(ax, means, stds, all_threads, modes)
    fig.tight_layout()
    save(fig, "plot_wall_time")

    # 2. Speedup
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_speedup(ax, means, all_threads, min_t, modes)
    fig.tight_layout()
    save(fig, "plot_speedup")

    # 3. Phase breakdown (one panel per mode)
    fig, ax_ph = plt.subplots(1, n_modes, figsize=(4 * n_modes + 1, 4), sharey=True)
    fig.suptitle("Phase time breakdown", fontsize=11, fontweight="bold")
    plot_phases(ax_ph if n_modes > 1 else [ax_ph], means, all_threads, modes)
    fig.tight_layout()
    save(fig, "plot_phases")

    # 4. Balance
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_balance(ax, parts, means, modes)
    fig.tight_layout()
    save(fig, "plot_balance")

    # 5. Peak RSS
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_rss(ax, means, all_threads, modes)
    fig.tight_layout()
    save(fig, "plot_rss")

    # ── Summary figure: 2-row GridSpec (top: wall+speedup | bottom: N×phases) ──
    fig = plt.figure(figsize=(4 * n_modes + 2, 9))
    fig.suptitle("tuna — partition strategy comparison", fontsize=13, fontweight="bold")

    # Bottom row: n_modes phase panels, each 2 columns wide → 2*n_modes total columns.
    # Top row: wall time (left half) + speedup (right half).
    gs = GridSpec(2, 2 * n_modes, figure=fig, hspace=0.42, wspace=0.35)
    ax_wall  = fig.add_subplot(gs[0, :n_modes])
    ax_speed = fig.add_subplot(gs[0, n_modes:])
    ax_ph    = [fig.add_subplot(gs[1, i*2:(i+1)*2]) for i in range(n_modes)]

    plot_wall_time(ax_wall,  means, stds, all_threads, modes)
    plot_speedup  (ax_speed, means, all_threads, min_t, modes)
    plot_phases   (ax_ph,    means, all_threads, modes)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, "plot_summary")

    print("\nDone.")


if __name__ == "__main__":
    main()
