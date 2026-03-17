#!/usr/bin/env bash
# bench_human.sh — tuna vs KMC on growing subsets of human genome files
#
# Usage: bash bench_human.sh <THREADS> <K> <KMC_RAM_GB>
# Example: bash bench_human.sh 4 31 200
#
# Input: FASTA.gz files listed in $FOF (multi-FASTA, one assembly per file)
# Sizes: N = 1, 2, 5, 10, 25, 50 files from the list
#
# Output: $RESULTS/bench.csv with per-run timing + unique k-mer counts
#         $RESULTS/tuna_n*.stderr  — tuna phase timings (phase1: / phase2:)
#         $RESULTS/kmc_n*.log      — KMC stdout (includes unique k-mer count)
#         $RESULTS/kmc_n*.time     — GNU time output for KMC count stage
#         $RESULTS/kmc_dump_n*.time— GNU time output for kmc_dump stage
#
# Correctness: unique k-mer count from tuna (wc -l) vs KMC stdout report.
# KMC dump is also run to get a full text output for any deeper comparison,
# but the dump wall time is reported separately since it is optional in practice.

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────

THREADS=${1:-4}
K=${2:-31}
KMC_RAM=${3:-200}       # GB — KMC -m flag

M=17                    # tuna minimizer length

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
FOF=/WORKS/vlevallois/data/dataset_genome_human/fof.list
WORK=/WORKS/vlevallois/test_tuna
RESULTS="$WORK/bench_human_$(date +%Y%m%d_%H%M%S)"

COUNTS=(1 2 5 10 25 50)

# ── Helpers ───────────────────────────────────────────────────────────────────

mkdir -p "$RESULTS"

# Convert GNU time "Elapsed (wall clock)" field (M:SS.ss or H:MM:SS.ss) → seconds
wall_to_s() {
    awk -F: '{if(NF==3) printf "%.2f", $1*3600+$2*60+$3;
              else       printf "%.2f", $1*60+$2}' <<< "$1"
}

# Extract max RSS in MB from GNU time -v output file
rss_mb() {
    local f="$1"
    local kb
    kb=$(grep "Maximum resident" "$f" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

# Extract wall time in seconds from GNU time -v output file
wall_from_file() {
    local f="$1"
    local t
    t=$(grep "Elapsed (wall clock)" "$f" | awk '{print $NF}')
    wall_to_s "$t"
}

# ── Header ────────────────────────────────────────────────────────────────────

CSV="$RESULTS/bench.csv"
echo "tool,n_files,wall_s,rss_mb,unique_kmers,phase1_s,phase2_s,kmc_dump_wall_s" > "$CSV"

echo "=== Human genome benchmark ==="
echo "    k=$K  m=$M  threads=$THREADS  kmc_ram=${KMC_RAM}GB"
echo "    fof: $FOF"
echo "    results: $RESULTS"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────

for N in "${COUNTS[@]}"; do

    echo "──── N = $N ──────────────────────────────────────────────────────────"

    SUBSET_FOF="$RESULTS/fof_${N}.list"
    head -n "$N" "$FOF" > "$SUBSET_FOF"

    # ── tuna ──────────────────────────────────────────────────────────────────

    TUNA_OUT="$RESULTS/tuna_n${N}.tsv"
    TUNA_WORK="$WORK/tuna_work_n${N}"
    mkdir -p "$TUNA_WORK"

    echo "  [tuna] running..."
    /usr/bin/time -v -o "$RESULTS/tuna_n${N}.time" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" \
        -w "$TUNA_WORK/" "@$SUBSET_FOF" "$TUNA_OUT" \
        2>"$RESULTS/tuna_n${N}.stderr" || {
            echo "  [tuna] FAILED — check $RESULTS/tuna_n${N}.stderr"
        }

    rm -rf "$TUNA_WORK"

    # Parse tuna results
    T_WALL=$(wall_from_file "$RESULTS/tuna_n${N}.time")
    T_RSS=$(rss_mb "$RESULTS/tuna_n${N}.time")
    T_KMERS=$([ -f "$TUNA_OUT" ] && wc -l < "$TUNA_OUT" || echo 0)
    T_P1=$(grep "^phase1:" "$RESULTS/tuna_n${N}.stderr" \
           | awk -F: '{gsub(/s/,"",$2); printf "%.2f", $2}')
    T_P2=$(grep "^phase2:" "$RESULTS/tuna_n${N}.stderr" \
           | awk -F: '{gsub(/s/,"",$2); printf "%.2f", $2}')

    echo "  [tuna] wall=${T_WALL}s  RSS=${T_RSS}MB  kmers=${T_KMERS}  p1=${T_P1}s  p2=${T_P2}s"
    echo "tuna,$N,$T_WALL,$T_RSS,$T_KMERS,$T_P1,$T_P2,na" >> "$CSV"

    # ── KMC ───────────────────────────────────────────────────────────────────

    KMC_DB="$WORK/kmc_db_n${N}"
    KMC_TMP="$WORK/kmc_tmp_n${N}"
    KMC_DUMP_OUT="$RESULTS/kmc_n${N}.tsv"
    mkdir -p "$KMC_TMP"

    echo "  [kmc]  counting..."
    /usr/bin/time -v -o "$RESULTS/kmc_n${N}.time" \
        kmc -k"$K" -m"$KMC_RAM" -ci1 -fm -t"$THREADS" \
        "@$SUBSET_FOF" "$KMC_DB" "$KMC_TMP" \
        > "$RESULTS/kmc_n${N}.log" 2>&1 || {
            echo "  [kmc]  FAILED — check $RESULTS/kmc_n${N}.log"
        }

    rm -rf "$KMC_TMP"

    echo "  [kmc]  dumping..."
    /usr/bin/time -v -o "$RESULTS/kmc_dump_n${N}.time" \
        kmc_dump -ci1 "$KMC_DB" "$KMC_DUMP_OUT" \
        >> "$RESULTS/kmc_n${N}.log" 2>&1 || {
            echo "  [kmc]  dump FAILED — check $RESULTS/kmc_n${N}.log"
        }

    rm -f "${KMC_DB}.kmc_pre" "${KMC_DB}.kmc_suf"

    # Parse KMC results
    K_WALL=$(wall_from_file "$RESULTS/kmc_n${N}.time")
    K_RSS=$(rss_mb "$RESULTS/kmc_n${N}.time")
    # KMC reports unique k-mers in its stdout log
    K_KMERS_LOG=$(grep -i "No. of unique counted k-mers" "$RESULTS/kmc_n${N}.log" \
                  | tail -1 | awk '{print $NF}' | tr -d ',')
    # Fallback: count dump lines
    if [ -z "$K_KMERS_LOG" ]; then
        K_KMERS_LOG=$([ -f "$KMC_DUMP_OUT" ] && wc -l < "$KMC_DUMP_OUT" || echo 0)
    fi
    K_DUMP_WALL=$(wall_from_file "$RESULTS/kmc_dump_n${N}.time")

    echo "  [kmc]  count_wall=${K_WALL}s  dump_wall=${K_DUMP_WALL}s  total=$(awk "BEGIN{printf \"%.2f\",$K_WALL+$K_DUMP_WALL}")s  RSS=${K_RSS}MB  kmers=${K_KMERS_LOG}"
    echo "kmc,$N,$K_WALL,$K_RSS,$K_KMERS_LOG,na,na,$K_DUMP_WALL" >> "$CSV"

    # ── Correctness ───────────────────────────────────────────────────────────

    if [ "$T_KMERS" -eq "$K_KMERS_LOG" ] 2>/dev/null; then
        echo "  [OK]   counts match: $T_KMERS unique k-mers"
    else
        DIFF=$(( T_KMERS - K_KMERS_LOG ))
        echo "  [DIFF] tuna=$T_KMERS  kmc=$K_KMERS_LOG  diff=$DIFF"
    fi

    echo ""
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=== Summary ==="
column -t -s, "$CSV"
echo ""
echo "Full results in: $RESULTS/"
