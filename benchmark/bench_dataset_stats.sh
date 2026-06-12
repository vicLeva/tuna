#!/usr/bin/env bash
# benchmark/bench_dataset_stats.sh — per-file sequence and k-mer statistics
#
# Processes each file individually (never accumulates across files).
# Uses a FIFO to pipe tuna output directly into awk so that the full k-mer
# TSV (up to ~100 GB for human) is never written to disk.
#
# Usage:
#   bash benchmark/bench_dataset_stats.sh [THREADS] [K]
#
# Output: $RESULTS/stats.csv, one row per file:
#   dataset, file_idx, filename,
#   n_seqs, total_bases, mean_len,   <- sequence stats (awk, no tool needed)
#   gc_pct,                          <- GC % of ACTG bases
#   distinct_k, total_k,             <- k-mer stats (tuna, k=K)
#   dup_rate, singleton_pct          <- total_k/distinct_k, fraction with count=1

set -euo pipefail

THREADS=${1:-8}
K=${2:-31}
M=21

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
WORK=/WORKS/vlevallois/expes_tuna
RESULTS="$WORK/dataset_stats_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS"

# ── Dataset registry ──────────────────────────────────────────────────────────
# Format: "name:source:fmt:n"
#   source  — fof path if n > 1, direct file path if n == 1
#   fmt     — fa (FASTA, plain or gz) | fq (FASTQ, plain or gz)
#   n       — number of files to process (head -n from fof, or 1)
DATASETS=(
    "ecoli:/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list:fa:100"
    "salmonella:/WORKS/vlevallois/data/dataset_pangenome_salmonella/fof.list:fa:100"
    "gut:/WORKS/vlevallois/data/dataset_metagenome_gut/fof.list:fa:100"
    "human:/WORKS/vlevallois/data/dataset_genome_human/fof.list:fa:10"
    "tara:/WORKS/vlevallois/data/dataset_metagenome_tara/fof.list:fq:10"
    "ecolir:/WORKS/vlevallois/data/data_sequencing/SRR2584863_1.fastq.gz:fq:1"
    "humanr:/WORKS/vlevallois/data/data_sequencing/SRR622461_1.fastq.gz:fq:1"
)

# ── CSV ───────────────────────────────────────────────────────────────────────
CSV="$RESULTS/stats.csv"
printf 'dataset,file_idx,filename,n_seqs,total_bases,mean_len,gc_pct,distinct_k,total_k,dup_rate,singleton_pct\n' \
    > "$CSV"
echo "Results: $CSV"
echo "k=$K  m=$M  threads=$THREADS"
echo ""

# ── Helpers ───────────────────────────────────────────────────────────────────

decompress() { [[ "$1" == *.gz ]] && gzip -dc "$1" || cat "$1"; }

# Sequence stats from FASTA (handles wrapped multi-line sequences + gzip).
# Outputs: n_seqs \t total_bases \t mean_len \t gc_pct
seq_stats_fa() {
    decompress "$1" | awk '
    BEGIN { n=0; total=0; gc=0; len=0; gci=0 }
    /^>/ {
        if (n > 0) { total += len; gc += gci }
        n++; len = 0; gci = 0; next
    }
    {
        len += length($0)
        gci += gsub(/[GCgc]/, "&")
    }
    END {
        if (n > 0) { total += len; gc += gci }
        mean = (n > 0) ? total / n : 0
        gcp  = (total > 0) ? gc * 100.0 / total : 0
        printf "%d\t%d\t%.1f\t%.2f\n", n, total, mean, gcp
    }'
}

# Sequence stats from FASTQ (every 4th line starting at 2, handles gzip).
# Outputs: n_reads \t total_bases \t mean_len \t gc_pct
seq_stats_fq() {
    decompress "$1" | awk '
    BEGIN { n=0; total=0; gc=0 }
    NR % 4 == 2 {
        l = length($0)
        total += l; n++
        gc += gsub(/[GCgc]/, "&")
    }
    END {
        mean = (n > 0) ? total / n : 0
        gcp  = (total > 0) ? gc * 100.0 / total : 0
        printf "%d\t%d\t%.1f\t%.2f\n", n, total, mean, gcp
    }'
}

# K-mer stats via tuna, streamed through a FIFO — no TSV written to disk.
# Outputs: distinct_k \t total_k \t dup_rate \t singleton_pct
kmer_stats() {
    local f="$1"
    local wdir="$RESULTS/tuna_work_$$_${RANDOM}"
    local fifo="$wdir/out.fifo"
    mkdir -p "$wdir"
    mkfifo "$fifo"

    # tuna writes the k-mer TSV to the FIFO; awk reads it concurrently.
    "$TUNA" -k "$K" -m "$M" -t "$THREADS" -hp -w "$wdir" "$f" "$fifo" \
        2>/dev/null &
    local tuna_pid=$!

    awk '
    BEGIN { d=0; total=0; s=0 }
    { d++; total += $2; if ($2 == 1) s++ }
    END {
        dup = (d > 0) ? total / d : 0
        sp  = (d > 0) ? s * 100.0 / d : 0
        printf "%d\t%d\t%.4f\t%.2f\n", d, total, dup, sp
    }' "$fifo"

    wait "$tuna_pid" || true
    rm -rf "$wdir"
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for ENTRY in "${DATASETS[@]}"; do
    IFS=: read -r DS SOURCE FMT N <<< "$ENTRY"

    echo "──── $DS  (fmt=$FMT, n=$N) ────────────────────────────"

    # Build file list
    if [ "$N" -eq 1 ]; then
        FILES=("$SOURCE")
    else
        if [ ! -f "$SOURCE" ]; then
            echo "  [SKIP] fof not found: $SOURCE"; continue
        fi
        mapfile -t FILES < <(head -n "$N" "$SOURCE")
    fi

    N_ACTUAL=${#FILES[@]}

    for IDX in "${!FILES[@]}"; do
        FILE="${FILES[$IDX]}"
        FNAME=$(basename "$FILE")
        FIDX=$(( IDX + 1 ))

        printf "  [%d/%d] %s\n" "$FIDX" "$N_ACTUAL" "$FNAME"

        # Sequence stats (no k-mer counting)
        printf "         seq_stats ..."
        if [ "$FMT" = "fa" ]; then
            SEQ=$(seq_stats_fa "$FILE")
        else
            SEQ=$(seq_stats_fq "$FILE")
        fi
        N_SEQS=$(printf '%s' "$SEQ" | cut -f1)
        TOTAL_BASES=$(printf '%s' "$SEQ" | cut -f2)
        MEAN_LEN=$(printf '%s' "$SEQ" | cut -f3)
        GC_PCT=$(printf '%s' "$SEQ" | cut -f4)
        printf " done (seqs=%s  bases=%s  mean=%.0fbp  gc=%s%%)\n" \
            "$N_SEQS" "$TOTAL_BASES" "$MEAN_LEN" "$GC_PCT"

        # K-mer stats via tuna
        printf "         kmer_stats ..."
        KS=$(kmer_stats "$FILE")
        DISTINCT=$(printf '%s' "$KS" | cut -f1)
        TOTAL_K=$(printf '%s' "$KS" | cut -f2)
        DUP=$(printf '%s' "$KS" | cut -f3)
        SING=$(printf '%s' "$KS" | cut -f4)
        printf " done (distinct=%s  total=%s  dup=%.2fx  singletons=%s%%)\n" \
            "$DISTINCT" "$TOTAL_K" "$DUP" "$SING"

        printf '%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$DS" "$FIDX" "$FNAME" \
            "$N_SEQS" "$TOTAL_BASES" "$MEAN_LEN" "$GC_PCT" \
            "$DISTINCT" "$TOTAL_K" "$DUP" "$SING" \
            >> "$CSV"
    done
    echo ""
done

echo "=== Done: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows) ==="
