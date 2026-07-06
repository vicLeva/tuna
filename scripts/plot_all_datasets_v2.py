#!/usr/bin/env python3
"""Per-file benchmark — tuna/KMC (old) + tuna_rob (rob_v2).

Reads:
  benchmark/results/bench_datasets_*/bench.csv  — most recent 3-tool CSV
  benchmark/results/rob_v2/bench_datasets.csv   — tuna_rob per-file CSV

Plots (dpi=500) in benchmark/results/rob_v2/:
  datasets_main.png  — 2 rows × 5 datasets (ecoli/salmonella/gut/human/tara)
                        row 0: wall-time boxplots (tuna dev · tuna_rob · KMC3)
                        row 1: per-file speedup tuna_rob / tuna_dev
"""

import csv, sys
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

RESULTS = Path(__file__).parent.parent / "benchmark" / "results"
ROB_CSV = RESULTS / "rob_v2" / "bench_datasets.csv"
OUT_DIR = RESULTS / "rob_v2"

C_TUNA    = "#2166ac"
C_ROB     = "#e6ab02"
C_KMC     = "#d6604d"
C_SPD_POS = "#59a14f"
C_SPD_NEG = C_KMC

DATASETS_ORDER = ["ecoli", "salmonella", "gut", "human", "tara"]


# ── loading ───────────────────────────────────────────────────────────────────

def _f(v):
    return float(v) if v not in ("", "na") else None

def _strip_s(v):
    """Strip trailing 's' from tuna stderr phase values (e.g. '0.051s' → 0.051)."""
    v = v.strip()
    if v.endswith("s"):
        v = v[:-1]
    return float(v) if v not in ("", "na", "") else None


def load_old(csv_path):
    """data[ds][fidx][tool] = {wall, p1, p2, rss}"""
    data = defaultdict(lambda: defaultdict(dict))
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, fidx, tool = row["dataset"], int(row["file_idx"]), row["tool"]
            if tool not in ("tuna", "kmc"):
                continue
            data[ds][fidx][tool] = dict(
                wall=_f(row["wall_s"]),
                p1=_f(row.get("phase1_s", "")),
                p2=_f(row.get("phase2_s", "")),
                rss=_f(row["rss_mb"]),
            )
    return data


def load_rob(csv_path):
    """data[ds][fidx] = {wall, p1, p2, rss}"""
    data = defaultdict(dict)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            ds, fidx = row["dataset"], int(row["file_idx"])
            data[ds][fidx] = dict(
                wall=_f(row["wall_s"]),
                p1=_strip_s(row.get("phase1_s", "")),
                p2=_strip_s(row.get("phase2_s", "")),
                rss=_f(row["rss_mb"]),
            )
    return data


# ── main plot ─────────────────────────────────────────────────────────────────

def plot_main(old, rob, out_path):
    datasets = [ds for ds in DATASETS_ORDER if ds in old and ds in rob]
    ND = len(datasets)

    fig, axes = plt.subplots(2, ND, figsize=(4 * ND, 9))
    if ND == 1:
        axes = axes.reshape(2, 1)

    fig.suptitle(
        "tuna (dev) vs tuna_rob vs KMC3 — per-file benchmark  (k=31, 8 threads)",
        fontsize=11, fontweight="bold")

    for col, ds in enumerate(datasets):
        # files present in all three sources
        common = sorted(
            set(old[ds].keys()) &
            set(k for k in old[ds] if "tuna" in old[ds][k]) &
            set(k for k in old[ds] if "kmc"  in old[ds][k]) &
            set(rob[ds].keys())
        )

        t_wall   = np.array([old[ds][i]["tuna"]["wall"] for i in common if old[ds][i].get("tuna",{}).get("wall") is not None])
        r_wall   = np.array([rob[ds][i]["wall"]         for i in common if rob[ds][i].get("wall")        is not None])
        k_wall   = np.array([old[ds][i]["kmc"]["wall"]  for i in common if old[ds][i].get("kmc", {}).get("wall")  is not None])

        # make lengths consistent (intersection of non-None across all three)
        valid = [i for i in common
                 if old[ds][i].get("tuna",{}).get("wall") is not None
                 and old[ds][i].get("kmc", {}).get("wall") is not None
                 and rob[ds][i].get("wall") is not None]
        t_wall = np.array([old[ds][i]["tuna"]["wall"] for i in valid])
        r_wall = np.array([rob[ds][i]["wall"]         for i in valid])
        k_wall = np.array([old[ds][i]["kmc"]["wall"]  for i in valid])
        n = len(valid)

        # ── Row 0: wall-time boxplots ──────────────────────────────────────
        ax = axes[0, col]
        bp = ax.boxplot([t_wall, r_wall, k_wall], patch_artist=True,
                        medianprops=dict(color="black", linewidth=2),
                        whiskerprops=dict(linewidth=1.2),
                        flierprops=dict(marker=".", markersize=4, alpha=0.5),
                        widths=0.5)
        for b, c in zip(bp["boxes"], [C_TUNA, C_ROB, C_KMC]):
            b.set_facecolor(c)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["tuna\n(dev)", "tuna_rob", "KMC3"], fontsize=8)
        ax.set_title(f"{ds}  (n={n})", fontweight="bold")
        ax.set_ylabel("Wall time (s)" if col == 0 else "")
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.grid(axis="y", alpha=0.3)
        ymax = ax.get_ylim()[1]
        for pos, vals, c in [(1, t_wall, C_TUNA), (2, r_wall, C_ROB), (3, k_wall, C_KMC)]:
            med = np.median(vals)
            ax.text(pos, ymax * 0.97, f"{med:.2f}s",
                    ha="center", va="top", fontsize=7.5, color=c, fontweight="bold")

        # ── Row 1: per-file speedup tuna_rob / tuna_dev ───────────────────
        ax = axes[1, col]
        speedup = t_wall / r_wall      # >1 means tuna_rob is faster
        order   = np.argsort(speedup)
        colors  = [C_SPD_POS if s >= 1 else C_SPD_NEG for s in speedup[order]]
        ax.bar(np.arange(n), speedup[order], color=colors, width=1.0)
        ax.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
        med = np.median(speedup)
        ax.axhline(med, color=C_ROB, linewidth=1.5, linestyle=":",
                   label=f"median {med:.2f}×")
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_xlabel("Files (sorted by speedup)")
        ax.set_ylabel("tuna_dev / tuna_rob  (>1 = rob faster)" if col == 0 else "")
        ax.set_title(f"tuna_rob speedup  med={med:.2f}×")
        ax.legend(fontsize=8, loc="lower right")
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=500, bbox_inches="tight")
    print(f"  saved  {out_path}")
    plt.close()


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    candidates = sorted(RESULTS.glob("bench_datasets_*/bench.csv"))
    if not candidates:
        sys.exit("error: no bench_datasets_*/bench.csv under benchmark/results/")
    old_csv = candidates[-1]

    print(f"Old CSV : {old_csv}")
    print(f"Rob CSV : {ROB_CSV}")
    print(f"Output  : {OUT_DIR}/datasets_main.png")

    old = load_old(old_csv)
    rob = load_rob(ROB_CSV)

    plot_main(old, rob, OUT_DIR / "datasets_main.png")


if __name__ == "__main__":
    main()
