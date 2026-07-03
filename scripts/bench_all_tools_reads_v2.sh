#!/usr/bin/env bash
# scripts/bench_all_tools_reads_v2.sh — Rob's tuna, reads datasets, thread sweep
#
# Mirrors bench_all_tools_reads.sh but tuna-only (Rob's fork).
# Datasets : G. gallus reads (SRX043656, 15 files) and
#            H. sapiens reads (ERR174324-ERR174341, 36 files)
#            Both passed as @fof (all files counted together).
# Threads  : configurable sweep (default: 1 2 4 6 8 10 16 22 28 32)
# m        : 21 for both (short reads; m=23 is for assembled genomes)
# n        : auto-tuned
#
# Output: $RESULTS_DIR/bench_reads.csv
#   dataset,n_files,threads,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers
# Aux files in $RESULTS_DIR/aux_reads/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"

RESULTS_DIR="/WORKS/vlevallois/expes_tuna/rob_v2"
WORKDIR="/WORKS/vlevallois/expes_tuna"

GALLUS_FOF="/WORKS/vlevallois/data/dataset_reads_gallus/fof.list"
HUMAN3_FOF="/WORKS/vlevallois/data/dataset_reads_human3/fof.list"

K=31
M=21
RAM_GB=256
THREADS_LIST=(1 2 4 6 8 10 16 22 28 32)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
[[ -z "$TUNA" ]] && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
[[ ! -f "$GALLUS_FOF" ]] && { echo "[error] not found: $GALLUS_FOF"; err=1; }
[[ ! -f "$HUMAN3_FOF" ]] && { echo "[error] not found: $HUMAN3_FOF"; err=1; }
[[ "$err" -eq 1 ]] && exit 1

GALLUS_N=$(wc -l < "$GALLUS_FOF")
HUMAN3_N=$(wc -l < "$HUMAN3_FOF")

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_reads"
CSV="$RESULTS_DIR/bench_reads.csv"
echo "dataset,n_files,threads,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers" > "$CSV"

echo "[bench] Rob's tuna — reads thread sweep"
echo "[bench] TUNA=$TUNA"
echo "[bench] k=$K  m=$M  ram=${RAM_GB}GB  threads=${THREADS_LIST[*]}"
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
# Run function
# =============================================================================
run_tuna() {
    local ds="$1" fof="$2" n_files="$3" threads="$4"
    local tag="${ds}_t${threads}"
    local tf="$RESULTS_DIR/aux_reads/${tag}.time"
    local se="$RESULTS_DIR/aux_reads/${tag}.stderr"
    local out="$WORKDIR/tuna_reads_${tag}.tsv"
    local work="$WORKDIR/tuna_reads_work_${tag}"
    mkdir -p "$work"

    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$threads" -ram "$RAM_GB" -hp \
        -w "$work/" "@$fof" "$out" \
        > /dev/null 2>"$se" \
    || { echo "  [FAIL] tuna $ds t=$threads"; rm -rf "$work" "$out"; return; }
    rm -rf "$work" "$out"

    local wall rss p1 p2 n_parts unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")
    n_parts=$(se_val "n_parts" "$se")
    unique=$(se_val "unique_kmers" "$se")

    printf "  %-8s n=%-3d t=%-3d m=%d  wall=%8.3fs  p1=%ss p2=%ss  RSS=%s MB  n=%s\n" \
        "$ds" "$n_files" "$threads" "$M" "$wall" "$p1" "$p2" "$rss" "$n_parts"
    echo "$ds,$n_files,$threads,$M,$wall,$rss,$p1,$p2,$n_parts,$unique" >> "$CSV"
}

# =============================================================================
# Main loop
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"
for threads in "${THREADS_LIST[@]}"; do
    echo "── t=$threads ──────────────────────────────────────────"
    run_tuna gallus  "$GALLUS_FOF"  "$GALLUS_N"  "$threads"
    run_tuna human3  "$HUMAN3_FOF"  "$HUMAN3_N"  "$threads"
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
