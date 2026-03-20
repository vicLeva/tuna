"""
tuna pipeline schematic — clean, no file names, no text overlaps.
Run: python3 pipeline_diagram.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# ── palette ───────────────────────────────────────────────────────────────────
C_IO   = "#4e79a7"
C_P1   = "#e07b28"
C_PART = "#3a9e96"
C_P2   = "#4a9e44"
C_OUT  = "#8f5fa8"
C_DARK = "#1a1a1a"
C_GRAY = "#555555"

fig, ax = plt.subplots(figsize=(17, 12))
ax.set_xlim(0, 17)
ax.set_ylim(0, 12)
ax.axis("off")
fig.patch.set_facecolor("white")
ax.set_facecolor("white")


# ── helpers ───────────────────────────────────────────────────────────────────

def rbox(ax, x, y, w, h, fc, title, lines=(),
         title_fs=9.5, line_fs=7.8, ec=None, zorder=3):
    ec = ec or C_DARK
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle="round,pad=0.05,rounding_size=0.18",
                       linewidth=1.4, edgecolor=ec, facecolor=fc, zorder=zorder)
    ax.add_patch(p)
    n = len(lines)
    cy = y + h / 2
    if n == 0:
        ax.text(x+w/2, cy, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color="white", zorder=zorder+1)
    else:
        gap = h / (n + 1.6)
        ty  = cy + gap * n / 2
        ax.text(x+w/2, ty, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color="white", zorder=zorder+1)
        for i, ln in enumerate(lines):
            ax.text(x+w/2, ty - (i+1)*gap, ln, ha="center", va="center",
                    fontsize=line_fs, color="white", alpha=0.90, zorder=zorder+1)

def sbox(ax, x, y, w, h, color, title, lines=(),
         title_fs=8.8, line_fs=7.4, zorder=4):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle="round,pad=0.04,rounding_size=0.12",
                       linewidth=1.1, edgecolor=color, facecolor="white",
                       alpha=0.95, zorder=zorder)
    ax.add_patch(p)
    n = len(lines)
    cy = y + h / 2
    if n == 0:
        ax.text(x+w/2, cy, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color=color, zorder=zorder+1)
    else:
        gap = h / (n + 1.8)
        ty  = cy + gap * n / 2
        ax.text(x+w/2, ty, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color=color, zorder=zorder+1)
        for i, ln in enumerate(lines):
            ax.text(x+w/2, ty - (i+1)*gap, ln, ha="center", va="center",
                    fontsize=line_fs, color=C_DARK, alpha=0.85, zorder=zorder+1)

def arrowh(ax, x0, x1, y, color=C_DARK, lw=1.8, label=""):
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color), zorder=8)
    if label:
        ax.text((x0+x1)/2, y+0.14, label, ha="center", va="bottom",
                fontsize=6.8, color=color, zorder=9)

def arrowv(ax, x, y0, y1, color=C_DARK, lw=1.8, label="", label_dx=0.14):
    ax.annotate("", xy=(x, y1), xytext=(x, y0),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color), zorder=8)
    if label:
        ax.text(x+label_dx, (y0+y1)/2, label, ha="left", va="center",
                fontsize=6.8, color=color, zorder=9)

def phase_line(ax, y, label, color):
    ax.plot([0.15, 16.85], [y, y], color=color, lw=0.9, ls="--", alpha=0.4, zorder=1)
    ax.text(0.2, y+0.08, label, fontsize=8.5, fontweight="bold",
            color=color, va="bottom", zorder=2)


# ══════════════════════════════════════════════════════════════════════════════
#  TITLE
# ══════════════════════════════════════════════════════════════════════════════
ax.text(8.5, 11.70, "tuna — Streaming k-mer Counter: Data Flow",
        ha="center", va="center", fontsize=16, fontweight="bold", color=C_DARK)
ax.text(8.5, 11.34, "FASTA / FASTQ  →  Phase 1: Partition  →  Phase 2: Count  →  k-mer TSV",
        ha="center", va="center", fontsize=9.5, color=C_GRAY)

# ══════════════════════════════════════════════════════════════════════════════
#  ROW 1 — INPUT  (y = 9.8 … 10.9)
# ══════════════════════════════════════════════════════════════════════════════
R1, RH1 = 9.80, 1.05

rbox(ax,  0.3,  R1, 2.5, RH1, C_IO, "Input Files",
     ("FASTA / FASTQ", "plain  or  .gz"))

rbox(ax,  3.5,  R1, 3.2, RH1, C_IO, "Sequence Reader",
     ("plain → mmap, zero-copy",
      ".gz  → 64 MB sliding window + decompress"))

rbox(ax,  7.5,  R1, 2.4, RH1, C_IO, "Strip non-ACTG",
     ("remove N, headers,", "newlines, whitespace"))

rbox(ax, 10.7,  R1, 2.6, RH1, C_IO, "ACTG-only chunk",
     ("contiguous DNA string,", "no ambiguity codes"))

arrowh(ax,  2.8,  3.5, R1+RH1/2, C_IO)
arrowh(ax,  6.7,  7.5, R1+RH1/2, C_IO)
arrowh(ax,  9.9, 10.7, R1+RH1/2, C_IO)

# threading note — placed to the right so it doesn't overlap the arrow
ax.text(16.8, R1 + RH1/2,
        "plain: threads\nsteal files\n\nsingle .gz:\n1 producer +\n(n−1) consumers",
        ha="right", va="center", fontsize=7.0, color=C_P1, style="italic",
        linespacing=1.4)

phase_line(ax, R1 - 0.22, "PHASE 1 — Partitioning", C_P1)

# ══════════════════════════════════════════════════════════════════════════════
#  ROW 2 — PHASE 1  (y = 7.2 … 9.55)
# ══════════════════════════════════════════════════════════════════════════════
R2, RH2 = 7.20, 2.35

# Background band
p = FancyBboxPatch((0.3, R2), 16.4, RH2,
                   boxstyle="round,pad=0.05,rounding_size=0.25",
                   linewidth=1.2, edgecolor=C_P1,
                   facecolor=C_P1, alpha=0.10, zorder=2)
ax.add_patch(p)

# Sub-box A: MinimizerWindow
sbox(ax, 0.55, R2+0.12, 3.5, RH2-0.24, C_P1,
     "Minimizer Window",
     ("ntHash rolling hash  (O(1) per base)",
      "canonical = fwd  XOR  rev",
      "two-stack sliding min  →  O(1)",
      "tracks minimizer position inline"))

# Sub-box B: Superkmer Extraction
sbox(ax, 4.35, R2+0.12, 3.8, RH2-0.24, C_P1,
     "Superkmer Extraction",
     ("flush superkmer on minimizer change",
      "or when length reaches 255 bases",
      "k−1 base overlap at boundaries",
      "record minimizer offset in superkmer"))

# Sub-box C: Superkmer Writer
sbox(ax, 8.45, R2+0.12, 4.9, RH2-0.24, C_P1,
     "Superkmer Writer  →  Binary Encoding",
     ("ASCII  →  2-bit packed  (4 bases / byte, MSB-first)",
      "header:  [ 1 B: length ]  [ 1 B: minimizer offset ]",
      "buffered per thread per partition",
      "flush under per-partition mutex"))

# arrows inside row 2
arrowh(ax, 4.05, 4.35, R2+RH2/2, C_P1, label="hash + pos")
arrowh(ax, 8.15, 8.45, R2+RH2/2, C_P1, label="superkmer")

# ACTG chunk drops into row 2
arrowv(ax, 12.0, R1, R2+RH2, C_P1, label="ACTG chunk")

# ══════════════════════════════════════════════════════════════════════════════
#  ROW 3 — PARTITION FILES  (y = 5.8 … 6.8)
# ══════════════════════════════════════════════════════════════════════════════
R3, RH3 = 5.80, 0.85

colors_part = [C_PART]*4 + [C_GRAY]
labels_part  = ["partition  0", "partition  1", "partition  2", "partition  3", "  …  "]
xs_part = [0.3, 3.5, 6.7, 9.9, 13.1]
ws_part = [2.9, 2.9, 2.9, 2.9, 1.6]

for lbl, xi, wi, c in zip(labels_part, xs_part, ws_part, colors_part):
    rbox(ax, xi, R3, wi, RH3, c, lbl, title_fs=9)

# arrow writer → files
arrowv(ax, 10.9, R2, R3+RH3, C_PART, label="flush()")

ax.text(8.5, R3 - 0.14,
        "n partitions  ·  auto-tuned:  next_pow2( total_input / 2 MB ),  clamped to  [16, fd_limit − 32]",
        ha="center", va="top", fontsize=7.5, color=C_GRAY)

phase_line(ax, R3 - 0.25, "PHASE 2 — Counting + Writing", C_P2)

# ══════════════════════════════════════════════════════════════════════════════
#  ROW 4 — PHASE 2  (y = 2.4 … 5.5)
# ══════════════════════════════════════════════════════════════════════════════
R4, RH4 = 2.40, 3.05

# Background band
p2 = FancyBboxPatch((0.3, R4), 16.4, RH4,
                    boxstyle="round,pad=0.05,rounding_size=0.25",
                    linewidth=1.2, edgecolor=C_P2,
                    facecolor=C_P2, alpha=0.10, zorder=2)
ax.add_patch(p2)

ax.text(8.5, R4+RH4-0.10,
        "one thread per partition  ·  no cross-partition locking",
        ha="center", va="top", fontsize=8.0, color=C_P2, style="italic", zorder=5)

# Sub-box A: Partition Reader
sbox(ax, 0.55, R4+0.12, 3.5, RH4-0.32, C_P2,
     "Partition Reader",
     ("mmap + MAP_POPULATE",
      "pages pre-faulted into RAM",
      "madvise SEQUENTIAL + WILLNEED",
      "next()  →  O(1), no syscalls",
      "parse 2-byte header inline"))

# Sub-box B: K-mer Counting
sbox(ax, 4.35, R4+0.12, 5.4, RH4-0.32, C_P2,
     "K-mer Counting",
     ("init: decode k bases  +  ntHash(l-mer at min_offset)  →  O(l)",
      "prefetch primary bucket  →  hide ~40 ns LLC miss",
      "upsert first k-mer of superkmer",
      "hot loop per remaining base:",
      "   2-bit unpack (bit-shift, no div)  →  advance window  →  upsert"))

# Sub-box C: Hash Table
sbox(ax, 10.05, R4+0.12, 3.5, RH4-0.32, C_P2,
     "Hash Table",
     ("flat bucket array  +  overflow table",
      "upsert:  key exists → count++",
      "         new key   → insert 1",
      "private per partition  (no lock)",
      "resize on load or overflow"))

# arrows inside row 4
arrowh(ax,  4.05,  4.35, R4+RH4/2, C_P2, label="packed bytes")
arrowh(ax,  9.75, 10.05, R4+RH4/2, C_P2, label="k-mer state")

# files → reader
arrowv(ax, 2.25, R3, R4+RH4, C_P2, label="mmap")

# ══════════════════════════════════════════════════════════════════════════════
#  ROW 5 — OUTPUT  (y = 0.5 … 1.8)
# ══════════════════════════════════════════════════════════════════════════════
R5, RH5 = 0.50, 1.25

rbox(ax, 4.2, R5, 8.0, RH5, C_OUT, "Output TSV",
     ("iterate all (k-mer, count) pairs  ·  filter:  ci  ≤  count  ≤  cx",
      "write  kmer\\tcount\\n  in 1 MB batches  ·  single output mutex"))

# hash table → output (curved arrow)
ax.annotate("", xy=(8.2, R5+RH5), xytext=(11.8, R4),
            arrowprops=dict(arrowstyle="-|>", lw=1.8, color=C_OUT,
                            connectionstyle="arc3,rad=-0.3"),
            zorder=8)

# ══════════════════════════════════════════════════════════════════════════════
#  LEGEND
# ══════════════════════════════════════════════════════════════════════════════
items = [
    mpatches.Patch(facecolor=C_IO,   edgecolor=C_DARK, label="I/O & parsing"),
    mpatches.Patch(facecolor=C_P1,   edgecolor=C_DARK, label="Phase 1: partitioning"),
    mpatches.Patch(facecolor=C_PART, edgecolor=C_DARK, label="Partition data"),
    mpatches.Patch(facecolor=C_P2,   edgecolor=C_DARK, label="Phase 2: counting"),
    mpatches.Patch(facecolor=C_OUT,  edgecolor=C_DARK, label="Output"),
]
ax.legend(handles=items, loc="lower right",
          bbox_to_anchor=(1.0, 0.0), fontsize=8.5,
          framealpha=0.95, edgecolor="#cccccc", ncol=5)

plt.tight_layout(pad=0.3)
plt.savefig("pipeline_diagram.png", dpi=180, bbox_inches="tight", facecolor="white")
plt.savefig("pipeline_diagram.pdf", bbox_inches="tight", facecolor="white")
print("Done.")
