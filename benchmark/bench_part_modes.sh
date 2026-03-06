#!/usr/bin/env bash
# ==============================================================================
# bench_part_modes.sh  —  tuna partition-strategy benchmark
#
# Compares the three partitioning modes tuna supports:
#   kmc      — KMC-style (default): norm-filtered canonical signature + table
#   kmtricks — kmtricks-style: pre-scan + load-balanced minimizer table  (-kmtricks)
#   hash     — simple hash: minimizer_hash % n_partitions               (-hash)
#
# For each mode × thread count × repetition the script records:
#   • total wall time and peak RSS (from GNU time -v)
#   • phase timings parsed from tuna's progress output:
#       phase01_s = scan (Phase 0, if any) + partition (Phase 1)
#       phase2_s  = count + write (Phase 2+3)
#   • partition load balance (per-file byte sizes via -kt flag)
#
# Outputs two CSV files in <out_dir>/:
#   runs.csv       — one row per (mode, threads, rep)
#   partitions.csv — one row per (run_id, partition_id)  [only rep=1 runs]
# and a per-run log:
#   <run_id>.log   — combined tuna stderr + GNU time report
#
# Requirements:
#   • tuna binary  (default: ../build/tuna/tuna)
#   • /usr/bin/time -v  (GNU time, NOT bash built-in time)
#   • grep -P  (PCRE-capable grep; standard on Linux)
#   • passwordless sudo recommended for page-cache drops
#
# Path overrides (export before calling, e.g. on a cluster):
#   TUNA_BIN=<path>   tuna binary      [default: ../build/tuna]
#   GNU_TIME=<path>   GNU time binary  [default: /usr/bin/time]
#
# Usage:
#   bash bench_part_modes.sh <fof_file> <out_dir> [KEY=VALUE ...]
#
# Key=Value options (can also be exported as env vars before the call):
#   K=<int>           k-mer length                   [default: 31]
#   M=<int>           minimizer length               [default: 17]
#   THREADS=<list>    quoted, space-separated list   [default: "1 2 4 8"]
#   PARTS=<int>       number of partitions           [default: 32]
#   REPS=<int>        repetitions per configuration  [default: 1]
#   KMC_SIG=<int>     KMC signature length (5-11)    [default: 9]
#
# Example:
#   bash bench_part_modes.sh ~/tmp/data/ecoli_fof_200.list results/ecoli200 \
#       K=31 THREADS="1 2 4 8" PARTS=32 REPS=3
# ==============================================================================
set -euo pipefail

# ── Paths (override via env vars if needed) ───────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TUNA="${TUNA_BIN:-$SCRIPT_DIR/../build/tuna}"
GNU_TIME="${GNU_TIME:-/usr/bin/time}"

# ── Parse positional arguments ────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    grep '^#' "$0" | sed 's/^# \?//'
    exit 1
fi

FOF="$1"
OUT_DIR="$2"
# If OUT_DIR has no path separator and is not absolute, place it under results/
[[ "$OUT_DIR" != */* && "$OUT_DIR" != /* ]] && OUT_DIR="$SCRIPT_DIR/results/$OUT_DIR"
shift 2

# ── Defaults (overridable via KEY=VALUE args or environment) ──────────────────
K="${K:-31}"
M="${M:-17}"
THREADS_LIST="${THREADS:-1 2 4 8}"
PARTS="${PARTS:-32}"
REPS="${REPS:-1}"
KMC_SIG="${KMC_SIG:-9}"

for arg in "$@"; do
    case "$arg" in
        K=*)        K="${arg#K=}" ;;
        M=*)        M="${arg#M=}" ;;
        THREADS=*)  THREADS_LIST="${arg#THREADS=}" ;;
        PARTS=*)    PARTS="${arg#PARTS=}" ;;
        REPS=*)     REPS="${arg#REPS=}" ;;
        KMC_SIG=*)  KMC_SIG="${arg#KMC_SIG=}" ;;
        *) echo "Unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ── Validate prerequisites ─────────────────────────────────────────────────────
[[ -x "$TUNA" ]]       || { echo "ERROR: tuna not found at $TUNA"          >&2; exit 1; }
[[ -f "$FOF"  ]]       || { echo "ERROR: file-of-files not found: $FOF"    >&2; exit 1; }
"$GNU_TIME" -v true 2>/dev/null \
    || { echo "ERROR: $GNU_TIME -v not available (need GNU time)" >&2; exit 1; }

# ── Output directories ────────────────────────────────────────────────────────
RUNS_CSV="$OUT_DIR/runs.csv"
PARTS_CSV="$OUT_DIR/partitions.csv"
SCRATCH="$OUT_DIR/.scratch"
mkdir -p "$OUT_DIR" "$SCRATCH"

# CSV headers
{
    printf 'run_id,mode,threads,partitions,rep'
    printf ',wall_s,peak_rss_kb'
    printf ',phase01_s,phase2_s,total_s'
    printf ',kmers_unique'
    printf ',p_mean_bytes,p_std_bytes,p_max_bytes,p_min_bytes,p_imbalance\n'
} > "$RUNS_CSV"

printf 'run_id,mode,threads,rep,partition_id,size_bytes\n' > "$PARTS_CSV"

# ── Page-cache helpers ────────────────────────────────────────────────────────
HAS_SUDO=false
sudo -n true 2>/dev/null && HAS_SUDO=true

drop_caches() {
    if $HAS_SUDO; then
        sync
        sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
    fi
}

# ── Parsing helpers ───────────────────────────────────────────────────────────

# Parse GNU time -v output from a combined log.
# Prints: wall_s  peak_rss_kb
parse_gnu_time() {
    local log="$1"
    local wall_raw wall_s rss
    wall_raw=$(grep "Elapsed (wall clock)" "$log" | awk '{print $NF}')
    wall_s=$(awk -F: '{
        if (NF == 3) printf "%.3f", ($1*3600 + $2*60 + $3)
        else         printf "%.3f", ($1*60  + $2)
    }' <<< "$wall_raw")
    rss=$(grep "Maximum resident" "$log" | awk '{print $NF}')
    echo "$wall_s $rss"
}

# Parse tuna's structured phase timings from a combined log.
# tuna always emits to stderr (regardless of -hp):
#   phase0: Xs   (only when pre-scan ran: kmc or kmtricks strategy)
#   phase1: Xs   (partitioning)
#   phase2: Xs   (count + write)
# Prints: phase01_s  phase2_s  total_s
# where phase01_s = phase0 + phase1 (total partitioning time incl. pre-scan).
parse_tuna_phases() {
    local log="$1"
    local p0 p1 p2 tot
    p0=$(grep  "^phase0:" "$log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
    p1=$(grep  "^phase1:" "$log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
    p2=$(grep  "^phase2:" "$log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
    tot=$(grep '\[done\] total:' "$log" 2>/dev/null \
        | grep -oP '[\d.]+(?=s)' || echo "0")
    p0="${p0:-0}"; p1="${p1:-0}"; p2="${p2:-0}"
    # Sum phase0 + phase1 for the combined partitioning column
    local p01; p01=$(awk "BEGIN{printf \"%.3f\", $p0 + $p1}")
    echo "$p01 $p2 $tot"
}

# Collect per-partition file sizes, append rows to PARTS_CSV, and print stats.
# Args: work_dir  n_parts  run_id  mode  threads  rep
# Prints: mean_bytes  std_bytes  max_bytes  min_bytes  imbalance
collect_partition_stats() {
    local wdir="$1" n="$2" run_id="$3" mode="$4" t="$5" rep="$6"
    local -a sizes
    for (( p = 0; p < n; ++p )); do
        local f="$wdir/${mode}_${p}.superkmers"
        local sz; sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
        sizes+=("$sz")
        # Only record individual partition rows on rep 1 (avoids huge CSV for multi-rep runs)
        if [[ "$rep" -eq 1 ]]; then
            printf '%s,%s,%s,%s,%d,%d\n' \
                "$run_id" "$mode" "$t" "$rep" "$p" "$sz" >> "$PARTS_CSV"
        fi
    done
    printf '%s\n' "${sizes[@]}" | awk '
        BEGIN { n=0; s=0; mx=0; mn=9999999999999 }
        { n++; s+=$1
          if ($1 > mx) mx=$1
          if ($1 < mn) mn=$1
          a[n]=$1 }
        END {
            mean = (n > 0) ? s/n : 0
            var  = 0
            for (i = 1; i <= n; i++) var += (a[i]-mean)^2
            std  = (n > 0) ? sqrt(var/n) : 0
            imb  = (mean > 0) ? mx/mean : 0
            printf "%.0f %.0f %.0f %.0f %.4f\n", mean, std, mx, mn, imb
        }'
}

# ── Mode definitions ──────────────────────────────────────────────────────────
MODES=(kmc kmtricks hash)

declare -A FLAGS
FLAGS[kmc]="-kmc -s $KMC_SIG"
FLAGS[kmtricks]="-kmtricks"
FLAGS[hash]="-hash"

# ── Summary header ────────────────────────────────────────────────────────────
n_threads_entries=$(echo "$THREADS_LIST" | wc -w)
total_runs=$(( ${#MODES[@]} * n_threads_entries * REPS ))

echo "============================================================"
echo "  bench_part_modes  —  tuna partition strategy benchmark"
printf "  k=%-3s  m=%-3s  partitions=%-4s  reps=%s\n" "$K" "$M" "$PARTS" "$REPS"
echo "  threads   : $THREADS_LIST"
echo "  kmc_sig   : $KMC_SIG"
echo "  input fof : $FOF  ($(wc -l < "$FOF") files)"
echo "  output    : $OUT_DIR"
echo "  binary    : $TUNA"
if $HAS_SUDO; then
    echo "  cache drops: ENABLED (fair, reproducible timings)"
else
    echo "  cache drops: DISABLED (no passwordless sudo — timings may be skewed)"
fi
echo "============================================================"
echo ""

# ── Main benchmark loop ───────────────────────────────────────────────────────
run_no=0

for mode in "${MODES[@]}"; do
    echo "── Mode: $mode ─────────────────────────────────────────────"

    for t in $THREADS_LIST; do
        for (( rep = 1; rep <= REPS; rep++ )); do
            (( run_no++ )) || true
            run_id="${mode}_t${t}_r${rep}"
            wdir="$SCRATCH/${run_id}"
            out_tsv="$SCRATCH/${run_id}.tsv"
            log="$OUT_DIR/${run_id}.log"

            mkdir -p "$wdir"

            printf "  [%2d/%2d] t=%2d  rep=%d  ... " \
                "$run_no" "$total_runs" "$t" "$rep"

            # Build command array — FLAGS may contain multiple words
            cmd=("$TUNA" -k "$K" -m "$M" -n "$PARTS" -t "$t" -kt -w "$wdir/")
            if [[ -n "${FLAGS[$mode]}" ]]; then
                read -ra extra_flags <<< "${FLAGS[$mode]}"
                cmd+=("${extra_flags[@]}")
            fi
            cmd+=("@$FOF" "$out_tsv")

            drop_caches

            # Run under GNU time; both tuna's stderr and the GNU time report
            # are written to the combined log (they use non-overlapping patterns).
            if ! "$GNU_TIME" -v "${cmd[@]}" 2>"$log"; then
                echo "FAILED (tuna exit non-zero — see $log)" >&2
                rm -rf "$wdir" "$out_tsv"
                continue
            fi

            # ── Collect metrics ──────────────────────────────────────────────
            read -r wall_s peak_rss_kb < <(parse_gnu_time "$log")
            read -r phase01_s phase2_s total_s < <(parse_tuna_phases "$log")
            kmers_unique=$(wc -l < "$out_tsv" 2>/dev/null || echo 0)

            # Partition balance (files persist because we used -kt)
            read -r p_mean p_std p_max p_min p_imbalance < <(
                collect_partition_stats "$wdir" "$PARTS" \
                    "$run_id" "$mode" "$t" "$rep"
            )

            # ── Append to runs.csv ───────────────────────────────────────────
            printf '%s,%s,%d,%d,%d,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s\n' \
                "$run_id" "$mode" "$t" "$PARTS" "$rep" \
                "$wall_s" "$peak_rss_kb" \
                "$phase01_s" "$phase2_s" "$total_s" \
                "$kmers_unique" \
                "$p_mean" "$p_std" "$p_max" "$p_min" "$p_imbalance" \
                >> "$RUNS_CSV"

            # ── Clean up scratch files for this run ──────────────────────────
            rm -rf "$wdir" "$out_tsv"

            printf "done  wall=%-7ss  RSS=%-8s kB  imbalance=%.2f\n" \
                "$wall_s" "$peak_rss_kb" "$p_imbalance"
        done
    done
    echo ""
done

# Remove empty scratch directory
rmdir "$SCRATCH" 2>/dev/null || true

# ── Final summary ─────────────────────────────────────────────────────────────
echo "============================================================"
echo "  Benchmark complete"
echo ""
echo "  CSV data:"
echo "    runs.csv       : $RUNS_CSV"
echo "    partitions.csv : $PARTS_CSV"
echo "  Per-run logs     : $OUT_DIR/<run_id>.log"
echo ""
echo "  Generate plots:"
echo "    python3 $(dirname "$0")/plot_part_modes.py $OUT_DIR"
echo "============================================================"
