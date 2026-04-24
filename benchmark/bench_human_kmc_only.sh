#!/usr/bin/env bash
# bench_human_kmc_only.sh — KMC count + dump on human genomes, matching bench_human_n_sweep.sh
# Outputs same CSV format; merge with existing tuna results locally.
#
# Usage: bash bench_human_kmc_only.sh [THREADS] [K] [M] [KMC_MEM_GB]

set -uo pipefail

THREADS=${1:-8}
K=${2:-31}
M=${3:-21}
KMC_MEM_GB=${4:-250}

REPS=3

KMC=${KMC:-kmc}
KMC_DUMP=${KMC_DUMP:-kmc_dump}
FOF=${FOF:-/WORKS/vlevallois/data/dataset_genome_human/fof.list}
WORK=${WORK:-/WORKS/vlevallois/test_tuna/human_nsweep}

N_VALUES=(1 2 3 5 7 10 15 20 30)

RESULTS="$WORK/kmc_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS"

TOTAL_FILES=$(wc -l < "$FOF")
echo "=== Human KMC-only (count+dump): k=$K  m=$M  t=$THREADS  mem=${KMC_MEM_GB}GB ==="
echo "    fof: $FOF  ($TOTAL_FILES files)"
echo "    results: $RESULTS"
echo ""

CSV="$RESULTS/bench_kmc.csv"
echo "tool,n_files,threads,rep,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers" > "$CSV"

wall_to_s() {
    awk -F: '{if(NF==3) printf "%.3f",$1*3600+$2*60+$3;
              else       printf "%.3f",$1*60+$2}' <<< "$1"
}
wall_from_file() { local t; t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}'); wall_to_s "$t"; }
rss_mb()         { local kb; kb=$(grep "Maximum resident" "$1" | awk '{print $NF}'); awk "BEGIN{printf \"%.0f\",$kb/1024}"; }

for N in "${N_VALUES[@]}"; do

    [ "$N" -gt "$TOTAL_FILES" ] && { echo "  [SKIP] N=$N > available ($TOTAL_FILES)"; continue; }

    SUBSET_FOF="$WORK/fof_n${N}.list"
    head -n "$N" "$FOF" > "$SUBSET_FOF"
    KMC_WORK="$WORK/kmc_only_work_n${N}"
    mkdir -p "$KMC_WORK"

    echo "── N=$N ──────────────────────────────────────────────────────────────────────"

    for rep in $(seq 1 $REPS); do
        stderr_f="$RESULTS/kmc_n${N}_r${rep}.stderr"
        time_f="$RESULTS/kmc_n${N}_r${rep}.timefile"

        /usr/bin/time -v -o "$time_f" bash -c "
            '$KMC' -k'$K' -t'$THREADS' -m'$KMC_MEM_GB' -fm -ci1 -hp \
                @'$SUBSET_FOF' '$KMC_WORK/out' '$KMC_WORK' &&
            '$KMC_DUMP' '$KMC_WORK/out' /dev/null
        " 2>"$stderr_f" || { echo "  [kmc FAIL N=$N rep=$rep]"; echo "kmc,$N,$THREADS,$rep,,,,,," >> "$CSV"; continue; }

        wall=$(wall_from_file "$time_f")
        rss=$(rss_mb "$time_f")
        s1=$(grep "^1st stage:" "$stderr_f" | awk '{print $3}' | tr -d 's')
        s2=$(grep "^2nd stage:" "$stderr_f" | awk '{print $3}' | tr -d 's')
        unique=$(grep "No. of unique counted k-mers" "$stderr_f" | awk '{print $NF}')

        echo "  [kmc] rep=$rep  wall=${wall}s  s1=${s1}s  s2=${s2}s  RSS=${rss}MB  unique=${unique}"
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
