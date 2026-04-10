#!/usr/bin/env bash
# bench_human_scale.sh — tuna vs KMC scaling with number of human genome files
#
# For each N in N_FILES_LIST, runs both tuna and KMC on the first N files
# of the human dataset, collecting timing, RSS, and (for tuna) per-partition
# table diagnostics.
#
# Usage: bash bench_human_scale.sh [THREADS] [K] [M] [KMC_RAM_GB]
# Example: bash bench_human_scale.sh 8 31 21 250
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
KMC_RAM=${4:-250}     # GB passed to KMC -m flag

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
KMC=kmc
KMC_DUMP=kmc_dump
FOF=/WORKS/vlevallois/data/dataset_genome_human/fof.list
WORK=/WORKS/vlevallois/test_tuna/human_scale
RESULTS="$WORK/results_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$RESULTS"

N_FILES_LIST=(1 2 3 5 10 15 20)

# ── Validate ──────────────────────────────────────────────────────────────────

[ -f "$FOF" ]    || { echo "ERROR: fof not found: $FOF"; exit 1; }
[ -x "$TUNA" ]   || { echo "ERROR: tuna binary not found: $TUNA"; exit 1; }
command -v "$KMC" >/dev/null || { echo "ERROR: kmc not in PATH"; exit 1; }

TOTAL_FILES=$(wc -l < "$FOF")
echo "=== Human genome scale experiment: tuna vs KMC ==="
echo "    k=$K  m=$M  threads=$THREADS  kmc_ram=${KMC_RAM}GB"
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

# ── run_kmc <tag> <fof_path> <n> ─────────────────────────────────────────────

run_kmc() {
    local tag="$1" fof="$2" n="$3"
    local kmc_db="$WORK/kmc_db_${tag}"
    local kmc_tmp="$WORK/kmc_tmp_${tag}"
    local kmc_log="$RESULTS/${tag}.kmc.log"
    local kmc_time_f="$RESULTS/${tag}.kmc.timefile"
    local dump_time_f="$RESULTS/${tag}.kmc_dump.timefile"

    mkdir -p "$kmc_tmp"

    # Phase 1: count k-mers.
    /usr/bin/time -v -o "$kmc_time_f" \
        "$KMC" -k"$K" -m"$KMC_RAM" -ci1 -fm -t"$THREADS" \
        "@${fof}" "$kmc_db" "$kmc_tmp" \
        >"$kmc_log" 2>&1 || {
            echo "  [kmc  FAIL] check $kmc_log"
            rm -rf "$kmc_tmp" "${kmc_db}.kmc_pre" "${kmc_db}.kmc_suf"
            return 1
        }
    rm -rf "$kmc_tmp"

    # Parse unique k-mer count from KMC log.
    local unique
    unique=$(grep -i "unique k-mers" "$kmc_log" | grep -v "below\|above" \
             | awk '{print $NF}' | head -1)
    [ -z "$unique" ] && unique=$(grep -i "No\. of unique" "$kmc_log" \
             | awk '{print $NF}' | head -1)
    [ -z "$unique" ] && unique=0

    # Phase 2: dump to /dev/null — measures time to scan and format all k-mers.
    /usr/bin/time -v -o "$dump_time_f" \
        "$KMC_DUMP" -ci1 "$kmc_db" /dev/null \
        >>"$kmc_log" 2>&1 || {
            echo "  [kmc_dump FAIL] check $kmc_log"
        }

    rm -f "${kmc_db}.kmc_pre" "${kmc_db}.kmc_suf"

    local p1 p2 wall rss
    p1=$(wall_from_file "$kmc_time_f")
    p2=$(wall_from_file "$dump_time_f")
    wall=$(awk "BEGIN{printf \"%.3f\", $p1 + $p2}")
    rss=$(rss_mb "$kmc_time_f")

    echo "  [kmc]   wall=${wall}s  count=${p1}s  dump=${p2}s  RSS=${rss}MB  unique=${unique}"

    echo "kmc,${n},${wall},${rss},${p1},${p2},${unique},\
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

    run_tuna "$TAG" "$TEMP_FOF" "$N" || true
    run_kmc  "$TAG" "$TEMP_FOF" "$N" || true

    rm -f "$TEMP_FOF"
    echo ""
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=== Done ==="
echo "Results: $RESULTS"
echo ""
column -t -s, "$CSV"
