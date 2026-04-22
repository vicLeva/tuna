#!/usr/bin/env bash
# bench_human_n_sweep.sh — tuna vs KMC, scaling with number of human genomes
#
# Usage: bash bench_human_n_sweep.sh [THREADS] [K] [M] [KMC_MEM_GB]
# Example: bash bench_human_n_sweep.sh 8 31 21 12
#
# Path overrides via env vars:
#   TUNA  KMC  FOF  WORK
#
# Output: $RESULTS/bench_n_sweep.csv
# Columns: tool, n_files, threads, rep, wall_s, rss_mb, phase1_s, phase2_s, unique_kmers

set -uo pipefail

THREADS=${1:-8}
K=${2:-31}
M=${3:-21}
KMC_MEM_GB=${4:-12}

REPS=3

# ── Paths ────────────────────────────────────────────────────────────────────
TUNA=${TUNA:-/WORKS/vlevallois/softs/tuna/build/tuna}
KMC=${KMC:-kmc}
FOF=${FOF:-/WORKS/vlevallois/data/dataset_genome_human/fof.list}
WORK=${WORK:-/WORKS/vlevallois/test_tuna/human_nsweep}

N_VALUES=(1 2 3 5 7 10 15 20 30)

RESULTS="$WORK/results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS"

TOTAL_FILES=$(wc -l < "$FOF")
echo "=== Human genome N-sweep: k=$K  m=$M  reps=$REPS  threads=$THREADS ==="
echo "    fof: $FOF  ($TOTAL_FILES files)"
echo "    results: $RESULTS"
echo ""

# ── CSV ───────────────────────────────────────────────────────────────────────
CSV="$RESULTS/bench_n_sweep.csv"
echo "tool,n_files,threads,rep,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers" > "$CSV"

# ── Helpers ───────────────────────────────────────────────────────────────────

wall_to_s() {
    awk -F: '{if(NF==3) printf "%.3f",$1*3600+$2*60+$3;
              else       printf "%.3f",$1*60+$2}' <<< "$1"
}

wall_from_file() {
    local t
    t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    wall_to_s "$t"
}

rss_mb() {
    local kb
    kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

parse_tuna() {
    grep "^${2}:" "$1" | awk -F': ' '{gsub(/s$/,"",$2); print $2}' | head -1
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for N in "${N_VALUES[@]}"; do

    if [ "$N" -gt "$TOTAL_FILES" ]; then
        echo "  [SKIP] N=$N > available ($TOTAL_FILES)"
        continue
    fi

    SUBSET_FOF="$WORK/fof_n${N}.list"
    head -n "$N" "$FOF" > "$SUBSET_FOF"

    echo "── N=$N ─────────────────────────────────────────────────────────────────────"

    # ── tuna ──────────────────────────────────────────────────────────────────
    TUNA_WORK="$WORK/tuna_work_n${N}"
    mkdir -p "$TUNA_WORK"

    for rep in $(seq 1 $REPS); do
        stderr_f="$RESULTS/tuna_n${N}_r${rep}.stderr"
        time_f="$RESULTS/tuna_n${N}_r${rep}.timefile"

        /usr/bin/time -v -o "$time_f" \
            "$TUNA" -k "$K" -m "$M" -t "$THREADS" -hp \
            -w "$TUNA_WORK/" "@$SUBSET_FOF" /dev/null \
            2>"$stderr_f" || {
                echo "  [tuna FAIL N=$N rep=$rep]"
                echo "tuna,$N,$THREADS,$rep,,,,,," >> "$CSV"
                continue
            }

        wall=$(wall_from_file "$time_f")
        rss=$(rss_mb "$time_f")
        p1=$(parse_tuna "$stderr_f" "phase1")
        p2=$(parse_tuna "$stderr_f" "phase2")
        unique=$(parse_tuna "$stderr_f" "unique_kmers")

        echo "  [tuna] rep=$rep  wall=${wall}s  p1=${p1}s  p2=${p2}s  RSS=${rss}MB  unique=${unique}"
        echo "tuna,$N,$THREADS,$rep,$wall,$rss,$p1,$p2,$unique" >> "$CSV"
    done
    rm -rf "$TUNA_WORK"

    # ── KMC ───────────────────────────────────────────────────────────────────
    KMC_WORK="$WORK/kmc_work_n${N}"
    mkdir -p "$KMC_WORK"

    for rep in $(seq 1 $REPS); do
        stderr_f="$RESULTS/kmc_n${N}_r${rep}.stderr"
        time_f="$RESULTS/kmc_n${N}_r${rep}.timefile"

        /usr/bin/time -v -o "$time_f" \
            "$KMC" -k"$K" -t"$THREADS" -m"$KMC_MEM_GB" -fm -ci1 -hp \
            "@$SUBSET_FOF" "$KMC_WORK/out" "$KMC_WORK" \
            2>"$stderr_f" || {
                echo "  [kmc FAIL N=$N rep=$rep]"
                echo "kmc,$N,$THREADS,$rep,,,,,," >> "$CSV"
                continue
            }

        wall=$(wall_from_file "$time_f")
        rss=$(rss_mb "$time_f")
        s1=$(grep "^1st stage:" "$stderr_f" | awk '{print $3}' | tr -d 's')
        s2=$(grep "^2nd stage:" "$stderr_f" | awk '{print $3}' | tr -d 's')
        unique=$(grep "No. of unique counted k-mers" "$stderr_f" | awk '{print $NF}')

        echo "  [kmc]  rep=$rep  wall=${wall}s  s1=${s1}s  s2=${s2}s  RSS=${rss}MB  unique=${unique}"
        echo "kmc,$N,$THREADS,$rep,$wall,$rss,$s1,$s2,$unique" >> "$CSV"
    done
    rm -rf "$KMC_WORK"

    rm -f "$SUBSET_FOF"
    echo ""

done

echo "=== Done ==="
echo "Results: $RESULTS"
echo ""
column -t -s, "$CSV"
