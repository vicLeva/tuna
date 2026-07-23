#!/usr/bin/env bash
# scripts/bench_big_data_dump.sh — large-scale reads benchmark: full dump
#
# Bis of bench_big_data.sh with real output written to disk instead of
# count-only, to measure the actual dump cost (vs. new_runs_v2's count-only
# baseline).
#
# One combined run per tool per dataset — all files in the dataset passed
# together (tuna/KMC via @fof, FastK as individual positional args), t=8
# fixed, m=21. Datasets: G. gallus reads (12 files) and H. sapiens reads
# (human3, 36 files).
#
# Full dump for all 3 tools — output files are deleted right after each
# run's timing is captured (not kept — a combined-dataset dump can be
# enormous, e.g. all 36 human3 files' k-mers in one TSV):
#   tuna  : writes its TSV (no -co). Dump cost is interleaved per-partition
#           into phase2 — compare phase2_s against new_runs_v2's count-only
#           phase2_s to estimate it.
#   KMC   : writes its binary db (no -w), then a separately-timed kmc_dump
#           step (dump_s) to ASCII. Compare 1st/2nd stage here against
#           new_runs_v2 to estimate KMC's own binary-write cost.
#   FastK : -t1 to produce the .ktab table, then Tabex -A to dump ASCII,
#           timed together as one combined wall time. Expect frequent
#           failures (Tabex crashes are common at this scale) — logged and
#           skipped, not fatal. Also still subject to the known
#           Scan_All_Input hang class (root-caused 2026-07-21, confirmed
#           independent of -t/-p) — kept under the same 6h timeout.
#
# Phase breakdown:
#   tuna : own "phase1:"/"phase2:" stderr lines
#   KMC  : "1st stage:"/"2nd stage:" stderr lines
#   FastK: no equivalent phase breakdown; left as "na"
#
# Output: $RESULTS_DIR/bench_big_data_dump.csv
#   dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_big_data/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"      # tuna binary (built from tuna's own repo, dev branch)
: "${KMC:=""}"       # kmc binary
: "${KMC_DUMP:=""}"  # kmc_dump binary
: "${FASTK:=""}"     # FastK binary
: "${TABEX:=""}"     # Tabex binary

EXPES_ROOT="/WORKS/vlevallois/expes_tuna"
RESULTS_DIR="$EXPES_ROOT/bench_big_data_dump_$(date +%Y%m%d_%H%M%S)"
WORKDIR="$EXPES_ROOT"

GALLUS_FOF="/WORKS/vlevallois/data/dataset_reads_gallus/fof.list"
HUMAN3_FOF="/WORKS/vlevallois/data/dataset_reads_human3/fof.list"

K=31
M=21
THREADS=8
RAM_GB=1536

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for var in TUNA KMC KMC_DUMP FASTK TABEX; do
    [[ -z "${!var}" ]]              && { echo "[error] $var is not set"; err=1; }
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] not executable: ${!var}"; err=1; }
done
[[ ! -f "$GALLUS_FOF" ]] && { echo "[error] not found: $GALLUS_FOF"; err=1; }
[[ ! -f "$HUMAN3_FOF" ]] && { echo "[error] not found: $HUMAN3_FOF"; err=1; }
[[ "$err" -eq 1 ]] && exit 1

GALLUS_N=$(wc -l < "$GALLUS_FOF")
HUMAN3_N=$(wc -l < "$HUMAN3_FOF")

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_big_data"
CSV="$RESULTS_DIR/bench_big_data_dump.csv"
echo "dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers" > "$CSV"

echo "[bench] big data benchmark — t=$THREADS, ram=${RAM_GB}GB, full dump"
echo "[bench] TUNA=$TUNA"
echo "[bench] KMC=$KMC"
echo "[bench] KMC_DUMP=$KMC_DUMP"
echo "[bench] FASTK=$FASTK"
echo "[bench] TABEX=$TABEX"
echo "[bench] k=$K  m=$M  threads=$THREADS"
echo "[bench] G. gallus  fof: $GALLUS_FOF ($GALLUS_N files)"
echo "[bench] H. sapiens fof: $HUMAN3_FOF ($HUMAN3_N files)"
echo "[bench] Results: $CSV"

# =============================================================================
# Helpers
# =============================================================================
wall_from_file() {
    local t; t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    awk -F: '{if(NF==3) printf "%.3f",$1*3600+$2*60+$3; else printf "%.3f",$1*60+$2}' <<< "$t"
}
rss_mb() {
    local kb; kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}
se_val() { grep "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2; gsub(/^ +| +s$/,"",v); print v}' || echo na; }

# =============================================================================
# Per-tool run functions
# Args: fof work timefile stderrfile [dumptimefile]
# =============================================================================

_run_tuna() {
    local fof="$1" work="$2" tf="$3" se="$4"
    # Full dump: real TSV output instead of /dev/null, no -co.
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" -hp \
        -w "$work/" "@$fof" "$work/out.tsv" \
        > /dev/null 2>"$se"
    local rc=$?
    rm -f "$work/out.tsv"
    return "$rc"
}

_run_kmc() {
    local fof="$1" work="$2" tf="$3" se="$4" dt="$5"
    # -fq: FASTQ input; @fof: list of files
    # Full dump: real binary db (no -w), then a separately-timed kmc_dump
    # step to get ASCII output. This cluster's KMC build prints its
    # 1st/2nd-stage summary to stdout, not stderr — capture both streams.
    mkdir -p "$work/tmp"
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k${K} -ci1 -cs4294967295 -fq -m${RAM_GB} -hp -t${THREADS} \
        "@$fof" "$work/out" "$work/tmp" \
        > "$se" 2>&1
    local rc=$?
    rm -rf "$work/tmp"
    if [[ "$rc" -ne 0 ]]; then return "$rc"; fi

    /usr/bin/time -v -o "$dt" \
        "$KMC_DUMP" "$work/out" "$work/out.tsv" \
        > /dev/null 2>>"$se"
    rc=$?
    rm -f "$work/out.tsv" "$work/out.kmc_pre" "$work/out.kmc_suf"
    return "$rc"
}

_run_fastk() {
    local fof="$1" work="$2" tf="$3" se="$4"
    # .fastq.gz is recognised natively by FastK; no symlinks needed
    # -t1: produce the .ktab table (count >= 1), then Tabex -A dumps ASCII.
    # Both timed together as one combined wall time. Expect frequent
    # failures here — Tabex crashes on a meaningful fraction of runs at
    # this scale.
    #
    # Known FastK bug (root-caused 2026-07-21): Scan_All_Input can hang
    # indefinitely on certain file-count/thread-count combinations
    # (confirmed N=threads+1 on a different dataset, and confirmed
    # independent of -t/-p; not a tuna issue). Wrapped in `timeout` so it
    # fails after 6h instead of stalling indefinitely.
    mkdir -p "$work/tmp"
    mapfile -t files < "$fof"
    timeout 21600 /usr/bin/time -v -o "$tf" bash -c "
        set -e
        \"$FASTK\" -k${K} -t1 -T${THREADS} -M${RAM_GB} \
            -N\"$work/out\" -P\"$work/tmp\" ${files[*]@Q} \
        && \"$TABEX\" -A \"$work/out\" > \"$work/out.tsv\"
    " > /dev/null 2>"$se"
    local rc=$?
    rm -rf "$work/tmp" "$work/out.tsv" "$work/out".*
    return "$rc"
}

# =============================================================================
# Dispatcher: run one tool on one dataset and record results
# =============================================================================
run_one() {
    local tool="$1" ds="$2" fof="$3" n_files="$4"
    local tag="${ds}_${tool}"
    local tf="$RESULTS_DIR/aux_big_data/${tag}.time"
    local se="$RESULTS_DIR/aux_big_data/${tag}.stderr"
    local dt="$RESULTS_DIR/aux_big_data/${tag}.dumptime"
    local work="$WORKDIR/big_data_dump_work_${tag}"
    mkdir -p "$work"

    local rc=0
    case "$tool" in
        tuna)  _run_tuna  "$fof" "$work" "$tf" "$se" ;;
        kmc)   _run_kmc   "$fof" "$work" "$tf" "$se" "$dt" ;;
        fastk) _run_fastk "$fof" "$work" "$tf" "$se" ;;
    esac
    rc=$?

    rm -rf "$work"
    if [[ "$rc" -eq 124 ]]; then
        echo "  [TIMEOUT] $tool $ds (known FastK N=threads+1-class hang, see memory)"
        echo "$ds,$n_files,$tool,timeout,,na,na,na,na" >> "$CSV"
        return
    fi
    [[ "$rc" -eq 0 ]] || { echo "  [FAIL] $tool $ds"; return; }

    local wall rss p1 p2 dump unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    case "$tool" in
        tuna)  p1=$(se_val "phase1"    "$se")
               p2=$(se_val "phase2"    "$se")
               dump=na
               unique=$(se_val "unique_kmers" "$se") ;;
        kmc)   p1=$(se_val "1st stage" "$se")
               p2=$(se_val "2nd stage" "$se")
               dump=$(wall_from_file "$dt")
               unique=$(grep "No. of unique k-mers" "$se" | awk -F: '{gsub(/^ +/,"",$2); print $2}') ;;
        fastk) p1=na; p2=na; dump=na; unique=na ;;
    esac

    printf "  %-6s %-8s n=%-3d  wall=%9ss  p1=%8ss  p2=%8ss  dump=%8ss  RSS=%8s MB  unique=%s\n" \
        "$tool" "$ds" "$n_files" "$wall" "$p1" "$p2" "$dump" "$rss" "$unique"
    echo "$ds,$n_files,$tool,$wall,$rss,$p1,$p2,$dump,$unique" >> "$CSV"
}

# =============================================================================
# Main loop: each tool, each dataset — one combined run per (tool, dataset)
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"
for tool in tuna kmc fastk; do
    echo "── $tool ──────────────────────────────────────────"
    run_one "$tool" gallus  "$GALLUS_FOF"  "$GALLUS_N"
    run_one "$tool" human3  "$HUMAN3_FOF"  "$HUMAN3_N"
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
