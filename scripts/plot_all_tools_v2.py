#!/usr/bin/env python3
"""Wall-time + RSS vs threads — tuna/KMC/FastK/KFC (old) + tuna_rob (rob_v2).

Reads:
  benchmark/results/results_all_tools.logs   — original 4-tool log
  benchmark/results/rob_v2/bench_tools.csv  — tuna_rob thread sweep CSV

Output PNGs (dpi=500) in benchmark/results/rob_v2/:
  tools_ecoli.png  tools_human.png  tools_combined.png
"""

import re, csv, sys, argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

RESULTS = Path(__file__).parent.parent / "benchmark" / "results"
OLD_LOG  = RESULTS / "results_all_tools.logs"
ROB_CSV  = RESULTS / "rob_v2" / "bench_tools.csv"
OUT_DIR  = RESULTS / "rob_v2"

THREADS = [1, 2, 4, 6, 8, 10, 16, 22, 28, 32]

TOOLS = ["tuna", "kmc", "fastk", "kfc", "tuna_rob"]
STYLE = {
    "tuna":     dict(color="#2166ac", marker="o",  label="tuna (dev)"),
    "kmc":      dict(color="#d6604d", marker="s",  label="KMC3"),
    "fastk":    dict(color="#4dac26", marker="^",  label="FastK"),
    "kfc":      dict(color="#8856a7", marker="D",  label="KFC"),
    "tuna_rob": dict(color="#e6ab02", marker="v",  label="tuna_rob", linewidth=2.2),
}

DATASET_META = {
    "ecoli": dict(title="E. coli  (k=31, 1 genome)", logtime=False, rss_gb=False),
    "human": dict(title="Human  (k=31, 1 genome)",   logtime=True,  rss_gb=True),
}


# ── parsing ───────────────────────────────────────────────────────────────────

def _wall_to_sec(s):
    parts = s.strip().split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def parse_log(path):
    hdr_re  = re.compile(r"=== TOOL=(\w+)\s+DS=(\w+)\s+FILE=\S+\s+T=(\d+)")
    wall_re = re.compile(r"Elapsed \(wall clock\) time.*?:\s+([\d:\.]+)")
    rss_re  = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
    exit_re = re.compile(r"Exit status:\s+(\d+)")

    text   = path.read_text(errors="replace")
    chunks = re.split(r"(?=^=== TOOL=)", text, flags=re.MULTILINE)

    records = []
    for chunk in chunks:
        m = hdr_re.search(chunk)
        if not m:
            continue
        tool, ds, threads = m.group(1).lower(), m.group(2).lower(), int(m.group(3))
        wm = wall_re.search(chunk)
        rm = rss_re.search(chunk)
        em = exit_re.search(chunk)
        records.append(dict(
            tool=tool, dataset=ds, threads=threads,
            wall_s=_wall_to_sec(wm.group(1)) if wm else None,
            rss_kb=int(rm.group(1))           if rm else None,
            ok=em is not None and int(em.group(1)) == 0,
        ))
    return records


def parse_rob_csv(path):
    records = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rss_mb = float(row["rss_mb"]) if row["rss_mb"] not in ("", "na") else None
            records.append(dict(
                tool="tuna_rob",
                dataset=row["dataset"],
                threads=int(row["threads"]),
                wall_s=float(row["wall_s"]) if row["wall_s"] not in ("", "na") else None,
                rss_kb=int(rss_mb * 1024) if rss_mb is not None else None,
                ok=True,
            ))
    return records


def lookup(records, tool, dataset, threads, field):
    hits = [r for r in records
            if r["tool"] == tool and r["dataset"] == dataset and r["threads"] == threads]
    if not hits:
        return None, None
    r = hits[0]
    return r.get(field), r.get("ok", True)


# ── formatters ────────────────────────────────────────────────────────────────

def _fmt_time(v, _):
    if v <= 0:   return "0 s"
    if v < 10:   return f"{v:.1f} s"
    return f"{v:.0f} s"

def _fmt_time_sparse(v, pos):
    import math
    if v <= 0: return ""
    exp = math.floor(math.log10(v))
    m = v / 10 ** exp
    if any(abs(m - x) < 0.05 for x in (1, 2, 5)):
        return _fmt_time(v, pos)
    return ""


# ── axis helpers ──────────────────────────────────────────────────────────────

def draw_time_ax(ax, records, dataset, logscale=False, add_legend=True):
    any_crash = False
    for tool in TOOLS:
        st = STYLE[tool]
        xs_ok, ys_ok   = [], []
        xs_bad, ys_bad = [], []
        for t in THREADS:
            val, ok = lookup(records, tool, dataset, t, "wall_s")
            if val is None:
                continue
            if ok:
                xs_ok.append(t);  ys_ok.append(val)
            else:
                xs_bad.append(t); ys_bad.append(val)
                any_crash = True
        lw = st.get("linewidth", 1.8)
        if xs_ok:
            ax.plot(xs_ok, ys_ok, color=st["color"], marker=st["marker"],
                    label=st["label"], linewidth=lw, markersize=6, zorder=3)
        if xs_bad:
            ax.plot(xs_bad, ys_bad, color=st["color"], marker="x",
                    linestyle="--", linewidth=1.2, markersize=9,
                    markeredgewidth=2.0, zorder=3,
                    label=f"{st['label']} (crashed)")
            if xs_ok and xs_bad[0] > xs_ok[-1]:
                ax.plot([xs_ok[-1], xs_bad[0]], [ys_ok[-1], ys_bad[0]],
                        color=st["color"], linestyle="--", linewidth=1.2, zorder=2)

    ax.set_xlabel("Threads", fontsize=9)
    ax.set_ylabel("Wall-clock time (log)" if logscale else "Wall-clock time", fontsize=9)
    ax.set_xticks(THREADS)
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    if logscale:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(
            mticker.LogLocator(base=10, subs=[1,2,3,4,5,6,7,8,9], numticks=20))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time_sparse))
        ax.yaxis.set_minor_locator(mticker.NullLocator())
        ax.set_ylim(bottom=60)
    else:
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    if any_crash:
        ax.annotate("✕ = crash", xy=(0.97, 0.03), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=7, color="grey", style="italic")
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85, loc="upper right")


def draw_rss_ax(ax, records, dataset, unit_gb=True, add_legend=False):
    unit_div   = 1024**2 if unit_gb else 1024
    unit_label = "GB" if unit_gb else "MB"
    for tool in TOOLS:
        st = STYLE[tool]
        xs, ys = [], []
        for t in THREADS:
            val, _ = lookup(records, tool, dataset, t, "rss_kb")
            if val is None:
                continue
            xs.append(t)
            ys.append(val / unit_div)
        lw = st.get("linewidth", 1.8)
        if xs:
            ax.plot(xs, ys, color=st["color"], marker=st["marker"],
                    label=st["label"], linewidth=lw, markersize=6, zorder=3)
    ax.set_xlabel("Threads", fontsize=9)
    ax.set_ylabel(f"Peak RSS ({unit_label})", fontsize=9)
    ax.set_xticks(THREADS)
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85)


# ── figures ───────────────────────────────────────────────────────────────────

def plot_dataset(records, dataset, meta, outprefix):
    fig, (ax_t, ax_r) = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    fig.suptitle(meta["title"], fontsize=11, fontweight="bold")
    draw_time_ax(ax_t, records, dataset, logscale=meta["logtime"], add_legend=True)
    draw_rss_ax(ax_r, records, dataset, unit_gb=meta["rss_gb"])
    out = f"{outprefix}_{dataset}.png"
    fig.savefig(out, dpi=500)
    plt.close(fig)
    print(f"  saved  {out}")


def plot_combined(records, outprefix):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting benchmark  (k=31, count + text dump)",
                 fontsize=11, fontweight="bold")

    for row, ds in enumerate(["ecoli", "human"]):
        meta = DATASET_META[ds]
        draw_time_ax(axes[row, 0], records, ds, logscale=False, add_legend=False)
        draw_rss_ax(axes[row, 1], records, ds, unit_gb=meta["rss_gb"], add_legend=False)
        rl = "E. coli" if ds == "ecoli" else "Human"
        axes[row, 0].set_ylabel(f"{rl}\nWall-clock time", fontsize=9)
        axes[row, 1].set_ylabel(f"Peak RSS ({'GB' if meta['rss_gb'] else 'MB'})", fontsize=9)

    # inset zoom on human wall-time
    ax_ht = axes[1, 0]
    axins = ax_ht.inset_axes([0.28, 0.24, 0.69, 0.56])
    for tool in ["tuna", "tuna_rob", "kmc", "fastk"]:
        st = STYLE[tool]
        xs, ys = [], []
        for t in THREADS:
            val, ok = lookup(records, tool, "human", t, "wall_s")
            if val is not None and ok:
                xs.append(t); ys.append(val)
        if xs:
            axins.plot(xs, ys, color=st["color"], marker=st["marker"],
                       linewidth=1.4, markersize=4, zorder=3)
    axins.set_xlim(THREADS[0] - 0.5, THREADS[-1] + 0.5)
    axins.set_ylim(60, 550)
    axins.set_xticks([1, 4, 8, 16, 22, 28])
    axins.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time))
    axins.tick_params(labelsize=6)
    axins.grid(axis="y", linestyle=":", alpha=0.35)
    axins.spines[["top", "right"]].set_visible(False)
    ax_ht.indicate_inset_zoom(axins, edgecolor="grey", alpha=0.55)

    seen, handles, labels = set(), [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in seen:
                handles.append(h); labels.append(l); seen.add(l)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.13, top=0.92)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.01),
               ncol=len(handles), fontsize=8.5, framealpha=0.85)

    out = f"{outprefix}_combined.png"
    fig.savefig(out, dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print(f"Old log : {OLD_LOG}")
    print(f"Rob CSV : {ROB_CSV}")
    print(f"Output  : {OUT_DIR}/tools_*.png")

    records = parse_log(OLD_LOG) + parse_rob_csv(ROB_CSV)
    print(f"Parsed {len(records)} run(s) total")

    outprefix = str(OUT_DIR / "tools")
    plot_combined(records, outprefix)


if __name__ == "__main__":
    main()
