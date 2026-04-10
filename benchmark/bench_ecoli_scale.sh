#!/usr/bin/env bash
# bench_ecoli_scale.sh — how does tuna performance scale with number of E. coli files?
#
# Hypothesis: the n_throughput floor in auto_tune_partitions grows O(N files)
# even when unique k-mer content doesn't (same species, high overlap).
# This causes over-partitioning and over-sized/under-used hash tables.
#
# Usage: bash bench_ecoli_scale.sh [THREADS] [K] [M]
# Example: bash bench_ecoli_scale.sh 8 31 21
#
# Output:
#   $RESULTS/bench_scale.csv        — one row per N_FILES run
#   $RESULTS/table_stats_N.csv      — per-partition table stats for each N
#   $RESULTS/<tag>.stderr           — full tuna stderr (timing + debug)
#
# All intermediate files go under $WORK — never /tmp.

set -euo pipefail

THREADS=${1:-8}
K=${2:-31}
M=${3:-21}

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
FOF=/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list
WORK=/WORKS/vlevallois/test_tuna/ecoli_scale
RESULTS="$WORK/results_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$RESULTS"

# Files to test — geometric progression covering 1..full dataset.
# Adjust or extend based on how many files exist in the fof.
N_FILES_LIST=(1 2 5 10 20 50 100 200 500 1000 2000 3600)

# ── Validate inputs ────────────────────────────────────────────────────────────

if [ ! -f "$FOF" ]; then
    echo "ERROR: fof not found: $FOF"
    exit 1
fi
if [ ! -x "$TUNA" ]; then
    echo "ERROR: tuna binary not found or not executable: $TUNA"
    exit 1
fi

TOTAL_FILES=$(wc -l < "$FOF")
echo "=== E. coli scale experiment ==="
echo "    k=$K  m=$M  threads=$THREADS"
echo "    fof: $FOF  ($TOTAL_FILES files total)"
echo "    results: $RESULTS"
echo ""

# ── CSV header ────────────────────────────────────────────────────────────────

CSV="$RESULTS/bench_scale.csv"
echo "n_files,n_parts,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers,\
dbg_load_mean,dbg_load_min,dbg_load_max,dbg_oversize_mean,dbg_unique_mean,\
dbg_n_resizes,dbg_parts_resized" > "$CSV"

# ── Helpers ───────────────────────────────────────────────────────────────────

wall_to_s() {
    awk -F: '{if(NF==3) printf "%.3f", $1*3600+$2*60+$3;
              else       printf "%.3f", $1*60+$2}' <<< "$1"
}

rss_mb() {
    local kb
    kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

parse_stderr() {
    local f="$1" key="$2"
    grep "^${key}:" "$f" | awk -F': ' '{gsub(/s$/,"",$2); print $2}' | head -1
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for N in "${N_FILES_LIST[@]}"; do

    # Cap at actual number of available files.
    if [ "$N" -gt "$TOTAL_FILES" ]; then
        echo "  [SKIP] N=$N > available files ($TOTAL_FILES)"
        continue
    fi

    TAG="ecoli_n$(printf '%05d' $N)"
    TUNA_WORK="$WORK/work_${TAG}"
    TEMP_FOF="$WORK/fof_${TAG}.list"

    # Build temp fof with first N files.
    mkdir -p "$TUNA_WORK"
    head -n "$N" "$FOF" > "$TEMP_FOF"

    echo "── N=$N files ──────────────────────────────────────────────────────────────"

    # Run tuna with -dbg to get per-partition stats, -kt to keep work dir for CSV.
    /usr/bin/time -v -o "$RESULTS/${TAG}.timefile" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" \
        -w "$TUNA_WORK/" -hp -dbg -kt \
        "@${TEMP_FOF}" /dev/null \
        2>"$RESULTS/${TAG}.stderr" || {
            echo "  [FAIL] tuna failed — check $RESULTS/${TAG}.stderr"
            rm -f "$TEMP_FOF"
            continue
        }

    # Parse timing from stderr.
    WALL=$(wall_to_s "$(grep 'Elapsed (wall clock)' "$RESULTS/${TAG}.timefile" | awk '{print $NF}')")
    RSS=$(rss_mb "$RESULTS/${TAG}.timefile")
    P1=$(parse_stderr "$RESULTS/${TAG}.stderr" "phase1")
    P2=$(parse_stderr "$RESULTS/${TAG}.stderr" "phase2")
    N_PARTS=$(parse_stderr "$RESULTS/${TAG}.stderr" "n_parts")
    UNIQUE=$(parse_stderr "$RESULTS/${TAG}.stderr" "unique_kmers")

    # Parse debug aggregate lines.
    DBG_LOAD_MEAN=$(parse_stderr "$RESULTS/${TAG}.stderr" "dbg_load_mean")
    DBG_LOAD_MIN=$(parse_stderr  "$RESULTS/${TAG}.stderr" "dbg_load_min")
    DBG_LOAD_MAX=$(parse_stderr  "$RESULTS/${TAG}.stderr" "dbg_load_max")
    DBG_OV_MEAN=$(parse_stderr   "$RESULTS/${TAG}.stderr" "dbg_oversize_mean")
    DBG_UQ_MEAN=$(parse_stderr   "$RESULTS/${TAG}.stderr" "dbg_unique_mean")
    DBG_RESIZES=$(parse_stderr   "$RESULTS/${TAG}.stderr" "dbg_n_resizes")
    DBG_PARTS_R=$(parse_stderr   "$RESULTS/${TAG}.stderr" "dbg_parts_resized")

    # Copy per-partition CSV from work dir.
    TABLE_CSV="$TUNA_WORK/debug_table_stats.csv"
    if [ -f "$TABLE_CSV" ]; then
        cp "$TABLE_CSV" "$RESULTS/table_stats_N${N}.csv"
    else
        echo "  [WARN] debug_table_stats.csv not found in $TUNA_WORK"
    fi

    # Cleanup.
    rm -rf "$TUNA_WORK" "$TEMP_FOF"

    echo "  n_parts=$N_PARTS  wall=${WALL}s  p1=${P1}s  p2=${P2}s  RSS=${RSS}MB"
    echo "  load: mean=${DBG_LOAD_MEAN}  min=${DBG_LOAD_MIN}  max=${DBG_LOAD_MAX}"
    echo "  oversize_mean=${DBG_OV_MEAN}  unique_mean=${DBG_UQ_MEAN}  resizes=${DBG_RESIZES}"

    echo "${N},${N_PARTS},${WALL},${RSS},${P1},${P2},${UNIQUE},\
${DBG_LOAD_MEAN},${DBG_LOAD_MIN},${DBG_LOAD_MAX},${DBG_OV_MEAN},${DBG_UQ_MEAN},\
${DBG_RESIZES},${DBG_PARTS_R}" >> "$CSV"

    echo ""
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=== Done ==="
echo "Results: $RESULTS"
echo ""
echo "bench_scale.csv:"
column -t -s, "$CSV"
