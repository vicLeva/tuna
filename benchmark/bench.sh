#!/usr/bin/env bash
# bench.sh — compare tuna vs KMC on a given file-of-files
#
# Usage:
#   bash bench.sh <fof_file> <run_name> [KEY=VALUE ...]
#
# Examples:
#   bash bench.sh /data/ecoli_fof_10.list ecoli10
#   bash bench.sh /data/human.list human K=31 M=17 THREADS=4 PARTS=128
#   bash bench.sh /data/tara.list tara FORMAT=fq K=31 THREADS=8 PARTS=128
#
# Key=Value options (can also be exported as env vars before the call):
#   K=<int>           k-mer length                       [default: 31]
#   M=<int>           minimizer length (must be < K)     [default: 17]
#   THREADS=<int>     worker threads                     [default: 1]
#   PARTS=<int>       number of tuna partitions          [default: 32]
#   STRATEGY=<str>    tuna partition strategy:           [default: hash]
#                       hash | kmtricks | kmc
#   FORMAT=<str>      input file format for KMC:        [default: fa]
#                       fa (FASTA) | fq (FASTQ)
#
# Path overrides (export before calling, e.g. on a cluster):
#   TUNA_BIN=<path>       tuna binary       [default: <script>/../build/tuna]
#   KMC_BIN=<path>        kmc binary        [default: <script>/../../kmc/bin/kmc]
#   KMC_TOOLS_BIN=<path>  kmc_tools binary  [default: <script>/../../kmc/bin/kmc_tools]
#   GNU_TIME=<path>       GNU time binary   [default: /usr/bin/time]
#   KMC_MEM=<int>         KMC RAM limit in GB            [default: 8]
#
# Scratch / tmp:
#   KMC scratch defaults to $TMPDIR if set (cluster local SSD), else $RUN_DIR/tmp.
#   Set TMPDIR to a fast local path on the cluster for best KMC performance.
#
# Page-cache fairness:
#   Each tool is given the same conditions by dropping the OS page cache
#   before its run (sync + echo 3 > /proc/sys/vm/drop_caches).  This requires
#   sudo.  If sudo is not available the cache is NOT dropped — results will
#   favour whichever tool runs second (warm cache), so comparisons across runs
#   may be unreliable.
set -euo pipefail

# ── Paths (override via env vars if needed) ───────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TUNA="${TUNA_BIN:-$SCRIPT_DIR/../build/tuna}"
KMC="${KMC_BIN:-$SCRIPT_DIR/../../kmc/bin/kmc}"
KMC_TOOLS="${KMC_TOOLS_BIN:-$SCRIPT_DIR/../../kmc/bin/kmc_tools}"
GNU_TIME="${GNU_TIME:-/usr/bin/time}"
TIMECMD="$GNU_TIME -v"
KMC_MEM="${KMC_MEM:-8}"

# ── Arguments ─────────────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    grep '^#' "$0" | sed 's/^# \?//'
    exit 1
fi

FOF="$1"
RUN_NAME="$2"
shift 2

# ── Parse KEY=VALUE overrides ─────────────────────────────────────────────────
for arg in "$@"; do
    case "$arg" in
        K=*)        K="${arg#K=}" ;;
        M=*)        M="${arg#M=}" ;;
        THREADS=*)  THREADS="${arg#THREADS=}" ;;
        PARTS=*)    PARTS="${arg#PARTS=}" ;;
        STRATEGY=*) STRATEGY="${arg#STRATEGY=}" ;;
        FORMAT=*)   FORMAT="${arg#FORMAT=}" ;;
        *) echo "tuna-bench: unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ── Parameters (env vars as fallback, then defaults) ──────────────────────────
K="${K:-31}"
M="${M:-17}"
THREADS="${THREADS:-1}"
PARTS="${PARTS:-32}"
STRATEGY="${STRATEGY:-hash}"
FORMAT="${FORMAT:-fa}"

# ── Validate ──────────────────────────────────────────────────────────────────
# Use command -v so bare names (e.g. KMC_BIN=kmc) are resolved via PATH.
require_bin() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $2 not found: $1" >&2; exit 1; }; }
require_bin "$TUNA"      "tuna"
require_bin "$KMC"       "kmc"
require_bin "$KMC_TOOLS" "kmc_tools"
[[ -f "$FOF" ]] || { echo "ERROR: file-of-files not found: $FOF" >&2; exit 1; }
"$GNU_TIME" -v true 2>/dev/null \
    || { echo "ERROR: $GNU_TIME -v not available (need GNU time, not bash builtin)" >&2; exit 1; }

case "$STRATEGY" in
    hash|kmtricks|kmc) ;;
    *) echo "ERROR: STRATEGY must be hash, kmtricks, or kmc (got: $STRATEGY)" >&2; exit 1 ;;
esac
case "$FORMAT" in
    fa|fq) ;;
    *) echo "ERROR: FORMAT must be fa or fq (got: $FORMAT)" >&2; exit 1 ;;
esac

# ── Output and scratch directories ────────────────────────────────────────────
RUN_DIR="$SCRIPT_DIR/results/$RUN_NAME"
# KMC scratch: prefer TMPDIR (cluster local SSD) over a sub-dir of RUN_DIR (may be NFS).
KMC_SCRATCH="${TMPDIR:-$RUN_DIR/tmp}/kmc_bench_$$"
mkdir -p "$RUN_DIR" "$KMC_SCRATCH"
trap 'rm -rf "$KMC_SCRATCH"' EXIT

# ── Build tuna strategy flag ──────────────────────────────────────────────────
case "$STRATEGY" in
    hash)     TUNA_STRATEGY_FLAGS=(-hash) ;;
    kmtricks) TUNA_STRATEGY_FLAGS=(-kmtricks) ;;
    kmc)      TUNA_STRATEGY_FLAGS=(-kmc) ;;
esac

# ── KMC format flag ───────────────────────────────────────────────────────────
# -fm = multi-FASTA file-of-files, -fq = FASTQ file-of-files
case "$FORMAT" in
    fa) KMC_FORMAT_FLAG="-fm" ;;
    fq) KMC_FORMAT_FLAG="-fq" ;;
esac

# ── Helpers ───────────────────────────────────────────────────────────────────

HAS_SUDO=false
sudo -n true 2>/dev/null && HAS_SUDO=true

drop_caches() {
    if $HAS_SUDO; then
        sync
        sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
    fi
}

# Convert h:mm:ss or m:ss.cc to fractional seconds
to_sec() {
    awk -F: '{
        if (NF == 3) printf "%.3f", ($1*3600 + $2*60 + $3)
        else         printf "%.3f", ($1*60  + $2)
    }'
}

# Extract a value from a GNU time -v log
gnu_time_val() {
    local log="$1" field="$2"
    case "$field" in
        wall) grep "Elapsed (wall clock)" "$log" | awk '{print $NF}' | to_sec ;;
        rss)  grep "Maximum resident"     "$log" | awk '{print $NF}' ;;
    esac
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "  Benchmark: $RUN_NAME"
echo "  k=$K   m=$M   threads=$THREADS   partitions=$PARTS"
echo "  tuna strategy: $STRATEGY"
echo "  input format:  $FORMAT"
echo "  input list:    $FOF  ($(wc -l < "$FOF") files)"
echo "  kmc RAM limit: ${KMC_MEM} GB"
if $HAS_SUDO; then
    echo "  cache drops: ENABLED (fair, reproducible timings)"
else
    echo "  cache drops: DISABLED (no passwordless sudo — timings may be skewed)"
fi
echo "============================================================"

# ── tuna ─────────────────────────────────────────────────────────────────────
echo ""
drop_caches
TUNA_CMD=("$TUNA"
    -k "$K" -m "$M"
    -n "$PARTS" -t "$THREADS"
    -hp
    "${TUNA_STRATEGY_FLAGS[@]}"
    "@$FOF"
    "$RUN_DIR/tuna_out.tsv")
echo "[tuna] ${TUNA_CMD[*]}"
$TIMECMD "${TUNA_CMD[@]}" 2>"$RUN_DIR/tuna_time.log"
echo "[tuna] done — $(wc -l < "$RUN_DIR/tuna_out.tsv") distinct k-mers"

# ── KMC (count) ───────────────────────────────────────────────────────────────
echo ""
drop_caches
KMC_COUNT_CMD=("$KMC"
    -k"$K"
    -m"$KMC_MEM"
    -ci1
    -t"$THREADS"
    "$KMC_FORMAT_FLAG"
    "@$FOF"
    "$KMC_SCRATCH/kmc_db"
    "$KMC_SCRATCH")
echo "[kmc]  ${KMC_COUNT_CMD[*]}"
$TIMECMD "${KMC_COUNT_CMD[@]}" 2>"$RUN_DIR/kmc_count_time.log"

# ── KMC (dump to TSV) ─────────────────────────────────────────────────────────
drop_caches
KMC_DUMP_CMD=("$KMC_TOOLS"
    transform
    "$KMC_SCRATCH/kmc_db"
    dump
    "$RUN_DIR/kmc_out.tsv")
echo "[kmc]  dump: ${KMC_DUMP_CMD[*]}"
$TIMECMD "${KMC_DUMP_CMD[@]}" 2>"$RUN_DIR/kmc_dump_time.log"
echo "[kmc]  done — $(wc -l < "$RUN_DIR/kmc_out.tsv") distinct k-mers"

# ── Summary ───────────────────────────────────────────────────────────────────
tuna_p0=$(grep "^phase0:" "$RUN_DIR/tuna_time.log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
tuna_p1=$(grep "^phase1:" "$RUN_DIR/tuna_time.log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
tuna_p2=$(grep "^phase2:" "$RUN_DIR/tuna_time.log" 2>/dev/null | awk '{gsub(/s$/,"",$2); print $2}' || true)
kmc_s1=$(grep  "^1st stage:" "$RUN_DIR/kmc_count_time.log" 2>/dev/null | awk '{gsub(/s$/,"",$3); print $3}' || true)
kmc_s2=$(grep  "^2nd stage:" "$RUN_DIR/kmc_count_time.log" 2>/dev/null | awk '{gsub(/s$/,"",$3); print $3}' || true)

tuna_wall=$(gnu_time_val "$RUN_DIR/tuna_time.log"           wall)
tuna_rss=$(gnu_time_val  "$RUN_DIR/tuna_time.log"           rss)
kmc_count_wall=$(gnu_time_val "$RUN_DIR/kmc_count_time.log" wall)
kmc_count_rss=$(gnu_time_val  "$RUN_DIR/kmc_count_time.log" rss)
kmc_dump_wall=$(gnu_time_val  "$RUN_DIR/kmc_dump_time.log"  wall)
kmc_dump_rss=$(gnu_time_val   "$RUN_DIR/kmc_dump_time.log"  rss)
kmc_total=$(awk "BEGIN{printf \"%.3f\", $kmc_count_wall + $kmc_dump_wall}")

tuna_kmers=$(wc -l < "$RUN_DIR/tuna_out.tsv")
kmc_kmers=$(wc -l  < "$RUN_DIR/kmc_out.tsv")

echo ""
echo "============================================================"
echo "  Results  (k=$K, m=$M, t=$THREADS, n=$PARTS, strategy=$STRATEGY)"
echo "============================================================"
printf "  k-mers written:  tuna=%-12s  kmc=%s\n" "$tuna_kmers" "$kmc_kmers"
echo ""
printf "  %-8s  %-26s  %9s  %10s\n" "tool" "phase" "time (s)" "peak RSS"
printf "  %-8s  %-26s  %9s  %10s\n" "--------" "--------------------------" "---------" "----------"
printf "  %-8s  %-26s  %9.3f  %7s kB\n" "tuna" "total (wall)"          "$tuna_wall"       "$tuna_rss"
if [[ -n "${tuna_p0:-}" ]]; then
printf "  %-8s  %-26s  %9s\n"           "tuna" "  phase0 (pre-scan)"    "$tuna_p0"
fi
printf "  %-8s  %-26s  %9s\n"           "tuna" "  phase1 (partition)"   "${tuna_p1:-N/A}"
printf "  %-8s  %-26s  %9s\n"           "tuna" "  phase2 (count+write)" "${tuna_p2:-N/A}"
echo ""
printf "  %-8s  %-26s  %9.3f  %7s kB\n" "kmc" "total (count+dump)"        "$kmc_total"       "$kmc_count_rss"
printf "  %-8s  %-26s  %9s\n"           "kmc" "  stage1 (partition)"       "${kmc_s1:-N/A}"
printf "  %-8s  %-26s  %9s\n"           "kmc" "  stage2 (count)"           "${kmc_s2:-N/A}"
printf "  %-8s  %-26s  %9.3f  %7s kB\n" "kmc" "  dump (binary→TSV)"        "$kmc_dump_wall"   "$kmc_dump_rss"
echo ""
echo "  Logs: $RUN_DIR/"
