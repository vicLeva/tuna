#!/usr/bin/env bash
# bench_human_scale.sh — tuna vs KMC scaling with number of human genome files
#
# For each N in N_FILES_LIST, runs both tuna and KMC on the first N files
# of the human dataset, collecting timing, RSS, and (for tuna) per-partition
# table diagnostics.
#
# Usage: bash bench_human_scale.sh [THREADS] [K] [M]
# Example: bash bench_human_scale.sh 8 31 21
# KMC results are hardcoded (results_20260410_145334). Re-run KMC manually if needed.
#
# Output (all under $RESULTS/):
#   bench_scale.csv              — one row per (tool, N)
#   table_stats_N<N>.csv         — tuna per-partition stats for each N
#   human_n<N>.{tuna,kmc}.stderr — raw tool output
#
# All intermediate files go under $WORK — never /tmp.

set -uo pipefail

THREADS=${1:-8}
K=${2:-31}
M=${3:-21}
TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
FOF=/WORKS/vlevallois/data/dataset_genome_human/fof.list
WORK=/WORKS/vlevallois/test_tuna/human_scale
RESULTS="$WORK/results_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$RESULTS"

N_FILES_LIST=(1 2 3 5 10 15 20)

# ── Validate ──────────────────────────────────────────────────────────────────

[ -f "$FOF" ]    || { echo "ERROR: fof not found: $FOF"; exit 1; }
[ -x "$TUNA" ]   || { echo "ERROR: tuna binary not found: $TUNA"; exit 1; }

TOTAL_FILES=$(wc -l < "$FOF")
echo "=== Human genome scale experiment: tuna vs KMC (KMC cached) ==="
echo "    k=$K  m=$M  threads=$THREADS"
echo "    fof: $FOF  ($TOTAL_FILES files total)"
echo "    results: $RESULTS"
echo ""

# ── CSV ───────────────────────────────────────────────────────────────────────

CSV="$RESULTS/bench_scale.csv"
echo "tool,n_files,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers,\
n_parts,dbg_load_mean,dbg_load_min,dbg_load_max,\
dbg_oversize_mean,dbg_unique_mean,dbg_n_resizes,dbg_parts_resized" > "$CSV"

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

wall_from_file() {
    local t
    t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    wall_to_s "$t"
}

parse_tuna() {   # parse_tuna <stderr_file> <key>   strips trailing 's'
    local f="$1" key="$2"
    grep "^${key}:" "$f" | awk -F': ' '{gsub(/s$/,"",$2); print $2}' | head -1
}

# ── run_tuna <tag> <fof_path> <n> ────────────────────────────────────────────

run_tuna() {
    local tag="$1" fof="$2" n="$3"
    local tuna_work="$WORK/tuna_work_${tag}"
    local stderr_f="$RESULTS/${tag}.tuna.stderr"
    local time_f="$RESULTS/${tag}.tuna.timefile"

    mkdir -p "$tuna_work"

    /usr/bin/time -v -o "$time_f" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" \
        -w "$tuna_work/" -hp -dbg -kt \
        "@${fof}" /dev/null \
        2>"$stderr_f" || {
            echo "  [tuna FAIL] check $stderr_f"
            rm -rf "$tuna_work"
            return 1
        }

    local wall rss p1 p2 n_parts unique
    wall=$(wall_from_file "$time_f")
    rss=$(rss_mb "$time_f")
    p1=$(parse_tuna "$stderr_f" "phase1")
    p2=$(parse_tuna "$stderr_f" "phase2")
    n_parts=$(parse_tuna "$stderr_f" "n_parts")
    unique=$(parse_tuna "$stderr_f" "unique_kmers")

    local load_mean load_min load_max ov_mean uq_mean n_resizes parts_resized
    load_mean=$(parse_tuna "$stderr_f" "dbg_load_mean")
    load_min=$(parse_tuna  "$stderr_f" "dbg_load_min")
    load_max=$(parse_tuna  "$stderr_f" "dbg_load_max")
    ov_mean=$(parse_tuna   "$stderr_f" "dbg_oversize_mean")
    uq_mean=$(parse_tuna   "$stderr_f" "dbg_unique_mean")
    n_resizes=$(parse_tuna "$stderr_f" "dbg_n_resizes")
    parts_resized=$(parse_tuna "$stderr_f" "dbg_parts_resized")

    # Copy per-partition CSV.
    local table_csv="$tuna_work/debug_table_stats.csv"
    [ -f "$table_csv" ] && cp "$table_csv" "$RESULTS/table_stats_N${n}.csv" \
        || echo "  [tuna] WARN: debug_table_stats.csv not found"

    rm -rf "$tuna_work"

    echo "  [tuna]  wall=${wall}s  p1=${p1}s  p2=${p2}s  RSS=${rss}MB  n_parts=${n_parts}  unique=${unique}"

    echo "tuna,${n},${wall},${rss},${p1},${p2},${unique},\
${n_parts},${load_mean},${load_min},${load_max},\
${ov_mean},${uq_mean},${n_resizes},${parts_resized}" >> "$CSV"
}

# ── KMC cached results (results_20260410_145334, k=31 m=21 t=8) ──────────────
# Hardcoded to avoid re-running KMC (each run costs ~10 min/file).
# Re-run KMC and update these if k, m, threads, or dataset changes.

declare -A KMC_WALL=( [1]=154.010 [2]=176.620 [3]=202.660 [5]=261.630 [10]=375.760 [15]=481.230 [20]=599.410 )
declare -A KMC_RSS=(  [1]=46741   [2]=91888   [3]=138214  [5]=222290  [10]=238471  [15]=238493  [20]=238529  )
declare -A KMC_P1=(   [1]=29.130  [2]=47.340  [3]=71.480  [5]=126.860 [10]=231.760 [15]=336.790 [20]=449.960 )
declare -A KMC_P2=(   [1]=124.880 [2]=129.280 [3]=131.180 [5]=134.770 [10]=144.000 [15]=144.440 [20]=149.450 )
declare -A KMC_UNQ=(  [1]=2496813618 [2]=2581788018 [3]=2625336481 [5]=2666269605 [10]=2724704476 [15]=2767924604 [20]=2796306069 )

emit_kmc_cached() {
    local n="$1"
    if [ -z "${KMC_WALL[$n]+x}" ]; then
        echo "  [kmc]   no cached result for N=$n — skipping"
        return
    fi
    echo "  [kmc]   wall=${KMC_WALL[$n]}s  count=${KMC_P1[$n]}s  dump=${KMC_P2[$n]}s  RSS=${KMC_RSS[$n]}MB  unique=${KMC_UNQ[$n]}  (cached)"
    echo "kmc,${n},${KMC_WALL[$n]},${KMC_RSS[$n]},${KMC_P1[$n]},${KMC_P2[$n]},${KMC_UNQ[$n]},\
na,na,na,na,na,na,na,na" >> "$CSV"
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for N in "${N_FILES_LIST[@]}"; do

    if [ "$N" -gt "$TOTAL_FILES" ]; then
        echo "  [SKIP] N=$N > available files ($TOTAL_FILES)"
        continue
    fi

    TAG="human_n$(printf '%03d' $N)"
    TEMP_FOF="$WORK/fof_${TAG}.list"
    head -n "$N" "$FOF" > "$TEMP_FOF"

    echo "── N=$N files ──────────────────────────────────────────────────────────────"

    run_tuna       "$TAG" "$TEMP_FOF" "$N" || true
    emit_kmc_cached "$N"

    rm -f "$TEMP_FOF"
    echo ""
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=== Done ==="
echo "Results: $RESULTS"
echo ""
column -t -s, "$CSV"
