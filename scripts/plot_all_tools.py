#!/usr/bin/env python3
"""Plot wall-time + peak RSS vs threads from bench_all_tools.sh's results.log.

Usage:
    python plot_all_tools.py <results_dir> [out.png]

<results_dir> is the directory produced by bench_all_tools.sh (e.g. after
scp'ing it back from the cluster) — must contain results.log.
"""

import re
import sys
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── per-tool visual identity ────────────────────────────────────────────────
TOOLS = ["tuna", "kmc", "fastk", "kfc"]
STYLE = {
    "tuna":  dict(color="#2166ac", marker="o", label="tuna"),
    "kmc":   dict(color="#d6604d", marker="s", label="KMC3"),
    "fastk": dict(color="#4dac26", marker="^", label="FastK"),
    "kfc":   dict(color="#8856a7", marker="D", label="KFC"),
}
THREADS = [1, 2, 4, 6, 8, 10, 16, 22, 28, 32]

DATASET_META = {
    "ecoli": dict(title="E. coli  (k=31, m=21, 1 genome)", logtime=False, rss_gb=False),
    "human": dict(title="Human  (k=31, m=21, 1 genome)",   logtime=True,  rss_gb=True),
}

# ── parsing ───────────────────────────────────────────────────────────────────


def _wall_to_sec(s):
    parts = s.strip().split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def parse_results(path):
    hdr_re = re.compile(r"=== TOOL=(\w+)\s+DS=(\w+)\s+FILE=\S+\s+T=(\d+)")
    wall_re = re.compile(r"Elapsed \(wall clock\) time.*?:\s+([\d:\.]+)")
    rss_re = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")
    exit_re = re.compile(r"Exit status:\s+(\d+)")

    text = path.read_text(errors="replace")
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
            rss_kb=int(rm.group(1)) if rm else None,
            ok=em is not None and int(em.group(1)) == 0,
        ))
    return records


def lookup(records, tool, dataset, threads, field):
    hits = [r for r in records
            if r["tool"] == tool and r["dataset"] == dataset and r["threads"] == threads]
    if not hits:
        return None, None
    r = hits[0]
    return r.get(field), r["ok"]


# ── formatters ────────────────────────────────────────────────────────────────


def _fmt_time(v, _):
    if v <= 0:
        return "0 s"
    if v < 10:
        return f"{v:.1f} s"
    return f"{v:.0f} s"


def _fmt_time_sparse(v, pos):
    if v <= 0:
        return ""
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
        xs_ok, ys_ok = [], []
        xs_bad, ys_bad = [], []
        for t in THREADS:
            val, ok = lookup(records, tool, dataset, t, "wall_s")
            if val is None:
                continue
            if ok:
                xs_ok.append(t); ys_ok.append(val)
            else:
                xs_bad.append(t); ys_bad.append(val)
                any_crash = True
        if xs_ok:
            ax.plot(xs_ok, ys_ok, color=st["color"], marker=st["marker"],
                    label=st["label"], linewidth=1.8, markersize=6, zorder=3)
        if xs_bad:
            ax.plot(xs_bad, ys_bad, color=st["color"], marker="x",
                    linestyle="--", linewidth=1.2, markersize=9,
                    markeredgewidth=2.0, zorder=3, label=f"{st['label']} (crashed)")
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
            mticker.LogLocator(base=10, subs=[1, 2, 3, 4, 5, 6, 7, 8, 9], numticks=20))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_time_sparse))
        ax.yaxis.set_minor_locator(mticker.NullLocator())
        ax.set_ylim(bottom=90)
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
    unit_div = 1024**2 if unit_gb else 1024
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
        if xs:
            ax.plot(xs, ys, color=st["color"], marker=st["marker"],
                    label=st["label"], linewidth=1.8, markersize=6, zorder=3)
    ax.set_xlabel("Threads", fontsize=9)
    ax.set_ylabel(f"Peak RSS ({unit_label})", fontsize=9)
    ax.set_xticks(THREADS)
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    if add_legend:
        ax.legend(fontsize=8, framealpha=0.85)


# ── combined figure ───────────────────────────────────────────────────────────


def plot_combined(records, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6))
    fig.suptitle("k-mer counting benchmark  (k=31, m=21, count + text dump)",
                 fontsize=11, fontweight="bold")

    for row, ds in enumerate(["ecoli", "human"]):
        meta = DATASET_META[ds]
        draw_time_ax(axes[row, 0], records, ds, logscale=meta["logtime"], add_legend=False)
        draw_rss_ax(axes[row, 1], records, ds, unit_gb=meta["rss_gb"], add_legend=False)
        rl = "E. coli" if ds == "ecoli" else "Human"
        axes[row, 0].set_ylabel(f"{rl}\nWall-clock time", fontsize=9)
        axes[row, 1].set_ylabel(f"Peak RSS ({'GB' if meta['rss_gb'] else 'MB'})", fontsize=9)

    # inset zoom on human wall-time
    ax_ht = axes[1, 0]
    axins = ax_ht.inset_axes([0.28, 0.24, 0.69, 0.56])
    for tool in TOOLS:
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
    axins.set_ylim(90, 550)
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

    fig.savefig(out_path, dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: plot_all_tools.py <results_dir> [out.png]")

    results_dir = Path(sys.argv[1])
    log_path = results_dir / "results.log"
    if not log_path.is_file():
        sys.exit(f"error: not found: {log_path}")

    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else results_dir / "tools_combined.png"

    records = parse_results(log_path)
    if not records:
        sys.exit("error: no run blocks found in the log file")
    print(f"Parsed {len(records)} run(s) from {log_path}")

    plot_combined(records, out_path)


if __name__ == "__main__":
    main()
