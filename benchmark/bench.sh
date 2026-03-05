#!/usr/bin/env bash
# bench.sh — compare tuna vs KMC on a given file-of-files
#
# Usage:
#   bash bench.sh <fof_file> <run_name> [K=<int>] [THREADS=<int>] [PARTS=<int>]
#
# Examples:
#   bash bench.sh /data/ecoli_fof_10.list ecoli10
#   bash bench.sh /data/human.list human K=31 THREADS=4 PARTS=128
#   K=31 THREADS=4 PARTS=128 bash bench.sh /data/human.list human
#
# Path overrides (export before calling, e.g. on a cluster):
#   TUNA_BIN=<path>       tuna binary       [default: ../build/tuna]
#   KMC_BIN=<path>        kmc binary        [default: ../../kmc/bin/kmc]
#   KMC_TOOLS_BIN=<path>  kmc_tools binary  [default: ../../kmc/bin/kmc_tools]
#   GNU_TIME=<path>       GNU time binary   [default: /usr/bin/time]
#
# Page-cache fairness:
#   Each tool is given the same conditions by dropping the OS page cache
#   before its run (sync + echo 3 > /proc/sys/vm/drop_caches).  This requires
#   sudo.  If sudo is not available the cache is NOT dropped — results will
#   favour whichever tool runs second (warm cache), so comparisons across runs
#   may be unreliable.  Run with sudo to get fair, reproducible timings.
set -euo pipefail

# ── Paths (override via env vars if needed) ───────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TUNA="${TUNA_BIN:-$SCRIPT_DIR/../build/tuna}"
KMC="${KMC_BIN:-$SCRIPT_DIR/../../kmc/bin/kmc}"
KMC_TOOLS="${KMC_TOOLS_BIN:-$SCRIPT_DIR/../../kmc/bin/kmc_tools}"
GNU_TIME="${GNU_TIME:-/usr/bin/time}"
TIMECMD="$GNU_TIME -v"

# ── Arguments ─────────────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <fof_file> <run_name> [K=<int>] [THREADS=<int>] [PARTS=<int>]" >&2
    exit 1
fi

FOF="$1"
RUN_NAME="$2"
shift 2

# Optional positional overrides: K=31 THREADS=4
for arg in "$@"; do
    case "$arg" in
        K=*)      K="${arg#K=}" ;;
        THREADS=*) THREADS="${arg#THREADS=}" ;;
        PARTS=*)  PARTS="${arg#PARTS=}" ;;
        *) echo "Unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ── Parameters (env vars as fallback) ─────────────────────────────────────────
K=${K:-31}
THREADS=${THREADS:-1}
PARTS=${PARTS:-32}

# ── Output directories ────────────────────────────────────────────────────────
RUN_DIR="$SCRIPT_DIR/results/$RUN_NAME"
WORK_DIR="$RUN_DIR/tmp"
mkdir -p "$RUN_DIR" "$WORK_DIR"

# ── Helpers ───────────────────────────────────────────────────────────────────

# Drop OS page cache if sudo is available; otherwise warn once.
HAS_SUDO=false
if sudo -n true 2>/dev/null; then
    HAS_SUDO=true
fi

drop_caches() {
    if $HAS_SUDO; then
        sync
        sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
    fi
}

# Parse GNU time -v output and print wall time + peak RSS
parse_time() {
    local log="$1" label="$2"
    local wall rss
    wall=$(grep "Elapsed (wall clock)" "$log" | awk '{print $NF}')
    rss=$(grep  "Maximum resident"     "$log" | awk '{print $NF}')
    printf "  %-30s  wall: %s   peak RSS: %s kB\n" "$label" "$wall" "$rss"
}

echo "============================================================"
echo "  Benchmark: $RUN_NAME"
echo "  k=$K   threads=$THREADS   partitions=$PARTS"
echo "  input list: $FOF"
if $HAS_SUDO; then
    echo "  cache drops: enabled (fair comparison)"
else
    echo "  cache drops: DISABLED (no passwordless sudo — timings may be skewed)"
fi
echo "============================================================"

# ── tuna ─────────────────────────────────────────────────────────────────────
echo ""
drop_caches
TUNA_CMD=("$TUNA" -k "$K" -n "$PARTS" -t "$THREADS" -hp @"$FOF" "$RUN_DIR/tuna_out.tsv")
echo "[tuna] cmd: ${TUNA_CMD[*]}"
$TIMECMD "${TUNA_CMD[@]}" \
    2>"$RUN_DIR/tuna_time.log"
echo "[tuna] done — $(wc -l < "$RUN_DIR/tuna_out.tsv") distinct k-mers written"

# ── KMC (count) ───────────────────────────────────────────────────────────────
echo ""
KMC_TMP="$WORK_DIR/kmc_tmp"
mkdir -p "$KMC_TMP"
drop_caches
KMC_COUNT_CMD=("$KMC" -k"$K" -ci1 -t"$THREADS" -fm @"$FOF" "$WORK_DIR/kmc_db" "$KMC_TMP")
echo "[kmc] cmd: ${KMC_COUNT_CMD[*]}"
$TIMECMD "${KMC_COUNT_CMD[@]}" \
    2>"$RUN_DIR/kmc_count_time.log"

# ── KMC (dump to TSV) ─────────────────────────────────────────────────────────
drop_caches
KMC_DUMP_CMD=("$KMC_TOOLS" transform "$WORK_DIR/kmc_db" dump "$RUN_DIR/kmc_out.tsv")
echo "[kmc] dump cmd: ${KMC_DUMP_CMD[*]}"
$TIMECMD "${KMC_DUMP_CMD[@]}" \
    2>"$RUN_DIR/kmc_dump_time.log"
echo "[kmc] done — $(wc -l < "$RUN_DIR/kmc_out.tsv") distinct k-mers written"

# ── Summary ───────────────────────────────────────────────────────────────────

# Convert h:mm:ss or m:ss.cc to fractional seconds
to_sec() {
    awk -F: '{
        if (NF == 3) print ($1*3600 + $2*60 + $3)
        else         print ($1*60  + $2)
    }'
}

# Extract a value from a GNU time log: field = "wall" or "rss"
gnu_time_val() {
    local log="$1" field="$2"
    case "$field" in
        wall) grep "Elapsed (wall clock)" "$log" | awk '{print $NF}' | to_sec ;;
        rss)  grep "Maximum resident"     "$log" | awk '{print $NF}' ;;
    esac
}

# Extract tuna and kmc phase times from stderr logs
tuna_p0=$(grep "^phase0:" "$RUN_DIR/tuna_time.log" | awk '{gsub(/s$/,"",$2); print $2}')
tuna_p1=$(grep "^phase1:" "$RUN_DIR/tuna_time.log" | awk '{gsub(/s$/,"",$2); print $2}')
tuna_p2=$(grep "^phase2:" "$RUN_DIR/tuna_time.log" | awk '{gsub(/s$/,"",$2); print $2}')
kmc_s1=$(grep  "^1st stage:" "$RUN_DIR/kmc_count_time.log" | awk '{gsub(/s$/,"",$3); print $3}')
kmc_s2=$(grep  "^2nd stage:" "$RUN_DIR/kmc_count_time.log" | awk '{gsub(/s$/,"",$3); print $3}')

tuna_wall=$(gnu_time_val "$RUN_DIR/tuna_time.log"      wall)
tuna_rss=$(gnu_time_val  "$RUN_DIR/tuna_time.log"      rss)
kmc_count_wall=$(gnu_time_val "$RUN_DIR/kmc_count_time.log" wall)
kmc_count_rss=$(gnu_time_val  "$RUN_DIR/kmc_count_time.log" rss)
kmc_dump_wall=$(gnu_time_val  "$RUN_DIR/kmc_dump_time.log"  wall)
kmc_dump_rss=$(gnu_time_val   "$RUN_DIR/kmc_dump_time.log"  rss)
kmc_total=$(echo "$kmc_count_wall + $kmc_dump_wall" | bc)

tuna_kmers=$(wc -l < "$RUN_DIR/tuna_out.tsv")
kmc_kmers=$(wc -l  < "$RUN_DIR/kmc_out.tsv")

echo ""
echo "============================================================"
echo "  Results  (k=$K, t=$THREADS, n=$PARTS)"
echo "============================================================"
printf "  k-mers written:  tuna=%-12s  kmc=%s\n" "$tuna_kmers" "$kmc_kmers"
echo ""
printf "  %-8s  %-24s  %9s  %10s\n" "tool" "phase" "time (s)" "peak RSS"
printf "  %-8s  %-24s  %9s  %10s\n" "--------" "------------------------" "---------" "----------"
printf "  %-8s  %-24s  %9.2f  %7s kB\n" "tuna" "total (wall)"        "$tuna_wall"       "$tuna_rss"
if [[ -n "${tuna_p0:-}" ]]; then
printf "  %-8s  %-24s  %9s\n"           "tuna" "phase0 (pre-scan)"    "$tuna_p0"
fi
printf "  %-8s  %-24s  %9s\n"           "tuna" "phase1 (partition)"   "${tuna_p1:-N/A}"
printf "  %-8s  %-24s  %9s\n"           "tuna" "phase2 (count+write)" "${tuna_p2:-N/A}"
echo ""
printf "  %-8s  %-24s  %9.2f  %7s kB\n" "kmc" "total (count+dump)"      "$kmc_total"       "$kmc_count_rss"
printf "  %-8s  %-24s  %9s\n"           "kmc" "stage1 (partition)"       "${kmc_s1:-N/A}"
printf "  %-8s  %-24s  %9s\n"           "kmc" "stage2 (radix-sort+count)" "${kmc_s2:-N/A}"
printf "  %-8s  %-24s  %9.2f  %7s kB\n" "kmc" "dump (binary→TSV)"        "$kmc_dump_wall"   "$kmc_dump_rss"
echo ""
echo "  Logs: $RUN_DIR/"

# ── Cleanup intermediate KMC files ───────────────────────────────────────────
rm -rf "$WORK_DIR"
