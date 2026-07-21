#!/usr/bin/env bash
# scripts/bench_big_data.sh — large-scale reads benchmark: tuna vs KMC vs FastK
#
# One combined run per tool per dataset — all files in the dataset passed
# together (tuna/KMC via @fof, FastK as individual positional args), t=8
# fixed, m=21. Datasets: G. gallus reads (12 files) and H. sapiens reads
# (human3, 36 files).
#
# Count-only for all 3 tools — no disk output, wall-time reflects counting
# alone (each tool's own native mechanism, verified against source/docs):
#   tuna  : -co (skip output writing after counting)
#   KMC   : -w (without output) — binary db never written
#   FastK : no -t/-p — DO_TABLE/DO_PROFILE stay 0, so Merge_Tables/
#           Merge_Profiles never run (core Sorting() still happens
#           unconditionally)
#
# Phase breakdown:
#   tuna : own "phase1:"/"phase2:" stderr lines
#   KMC  : "1st stage:"/"2nd stage:" stderr lines — printed unconditionally
#          by kmc_CLI/kmc.cpp's print_summary() (undocumented in -h usage
#          text, but always emitted; verified against KMC's own CLI source
#          and empirically against the local binary)
#   FastK: no equivalent phase breakdown found; left as "na"
#
# Output: $RESULTS_DIR/bench_big_data.csv
#   dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_big_data/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"      # tuna binary (built from tuna's own repo, dev branch)
: "${KMC:=""}"       # kmc binary
: "${FASTK:=""}"     # FastK binary

EXPES_ROOT="/WORKS/vlevallois/expes_tuna"
RESULTS_DIR="$EXPES_ROOT/bench_big_data_$(date +%Y%m%d_%H%M%S)"
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
for var in TUNA KMC FASTK; do
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
CSV="$RESULTS_DIR/bench_big_data.csv"
echo "dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers" > "$CSV"

echo "[bench] big data benchmark — t=$THREADS, ram=${RAM_GB}GB, count-only"
echo "[bench] TUNA=$TUNA"
echo "[bench] KMC=$KMC"
echo "[bench] FASTK=$FASTK"
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
# Args: fof work timefile stderrfile
# =============================================================================

_run_tuna() {
    local fof="$1" work="$2" tf="$3" se="$4"
    # -co: skip output writing after counting
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" -hp -co \
        -w "$work/" "@$fof" /dev/null \
        > /dev/null 2>"$se"
}

_run_kmc() {
    local fof="$1" work="$2" tf="$3" se="$4"
    # -fq: FASTQ input; @fof: list of files
    # -w: without output — binary db is never written (true count-only)
    # Note: this cluster's KMC build prints its 1st/2nd-stage summary to
    # stdout, not stderr — capture both streams into $se so se_val can find it.
    mkdir -p "$work/tmp"
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k${K} -ci1 -cs4294967295 -fq -m${RAM_GB} -hp -t${THREADS} -w \
        "@$fof" "$work/out" "$work/tmp" \
        > "$se" 2>&1
    rm -rf "$work/tmp"
}

_run_fastk() {
    local fof="$1" work="$2" tf="$3" se="$4"
    # .fastq.gz is recognised natively by FastK; no symlinks needed
    # No -t/-p: DO_TABLE/DO_PROFILE stay 0, so no .ktab ever gets written —
    # true count-only (core Sorting() still happens unconditionally either way)
    mkdir -p "$work/tmp"
    mapfile -t files < "$fof"
    /usr/bin/time -v -o "$tf" \
        "$FASTK" -k${K} -T${THREADS} -M${RAM_GB} \
        -N"$work/out" -P"$work/tmp" \
        "${files[@]}" \
        > /dev/null 2>"$se"
    rm -rf "$work/tmp"
}

# =============================================================================
# Dispatcher: run one tool on one dataset and record results
# =============================================================================
run_one() {
    local tool="$1" ds="$2" fof="$3" n_files="$4"
    local tag="${ds}_${tool}"
    local tf="$RESULTS_DIR/aux_big_data/${tag}.time"
    local se="$RESULTS_DIR/aux_big_data/${tag}.stderr"
    local work="$WORKDIR/big_data_work_${tag}"
    mkdir -p "$work"

    local ok=true
    case "$tool" in
        tuna)  _run_tuna  "$fof" "$work" "$tf" "$se" ;;
        kmc)   _run_kmc   "$fof" "$work" "$tf" "$se" ;;
        fastk) _run_fastk "$fof" "$work" "$tf" "$se" ;;
    esac || ok=false

    rm -rf "$work"
    $ok || { echo "  [FAIL] $tool $ds"; return; }

    local wall rss p1 p2 unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    case "$tool" in
        tuna)  p1=$(se_val "phase1"    "$se")
               p2=$(se_val "phase2"    "$se")
               unique=$(se_val "unique_kmers" "$se") ;;
        kmc)   p1=$(se_val "1st stage" "$se")
               p2=$(se_val "2nd stage" "$se")
               unique=$(grep "No. of unique k-mers" "$se" | awk -F: '{gsub(/^ +/,"",$2); print $2}') ;;
        fastk) p1=na; p2=na; unique=na ;;
    esac

    printf "  %-6s %-8s n=%-3d  wall=%9ss  p1=%8ss  p2=%8ss  RSS=%8s MB  unique=%s\n" \
        "$tool" "$ds" "$n_files" "$wall" "$p1" "$p2" "$rss" "$unique"
    echo "$ds,$n_files,$tool,$wall,$rss,$p1,$p2,$unique" >> "$CSV"
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
