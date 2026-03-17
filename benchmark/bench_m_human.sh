#!/usr/bin/env bash
# bench_m_human.sh — impact of minimizer length m on a human genome dataset
#
# Usage: bash bench_m_human.sh <N_FILES> <THREADS> <N_PARTS>
# Example: bash bench_m_human.sh 1 4 4096
#
# Path overrides via env vars (defaults = cluster paths):
#   TUNA=/path/to/tuna  FOF=/path/to/fof.list  WORK=/path/to/workdir
# Local example:
#   TUNA=~/softs/tuna/build/tuna FOF=~/data/human.gz WORK=~/tuna_bench_tmp \
#     bash bench_m_human.sh 1 4 4096
#
# Runs tuna with m in {17,19,21,23,25,27} for k=31.
# Reports phase1/phase2 times, overflow stats, and resize stats (-dbg) per m.

set -euo pipefail

N_FILES=${1:-1}
THREADS=${2:-4}
N_PARTS=${3:-4096}

K=31
M_VALUES=(15 17 19 21 23 25 27)

TUNA=${TUNA:-/WORKS/vlevallois/softs/tuna/build/tuna}
FOF=${FOF:-/WORKS/vlevallois/data/dataset_genome_human/fof.list}
WORK=${WORK:-/WORKS/vlevallois/test_tuna}
RESULTS="$WORK/bench_m_human_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$RESULTS"

# Build subset fof: if FOF is a plain file (not a list), use it directly.
SUBSET_FOF="$RESULTS/fof_${N_FILES}.list"
if [[ "$FOF" == *.list || "$FOF" == *.txt ]]; then
    head -n "$N_FILES" "$FOF" > "$SUBSET_FOF"
else
    # Single file passed directly — write it as a 1-entry list.
    echo "$FOF" > "$SUBSET_FOF"
fi

CSV="$RESULTS/bench_m.csv"
echo "m,phase1_s,phase2_s,total_s,overflow_pct,overflow_kmers,partitions,n_resizes,resize_s,n_ov_resize,n_load_resize" > "$CSV"

echo "=== m-value benchmark: k=$K  n_files=$N_FILES  threads=$THREADS  n_parts=$N_PARTS ==="
echo "Results: $RESULTS"
echo ""
printf "%-4s  %-9s  %-9s  %-9s  %-10s  %-14s  %-9s  %-9s  %-10s\n" \
    "m" "phase1_s" "phase2_s" "total_s" "overflow_%" "overflow_kmers" "n_resizes" "resize_s" "ov/load_trig"
echo "────  ─────────  ─────────  ─────────  ──────────  ──────────────  ─────────  ─────────  ──────────"

for M in "${M_VALUES[@]}"; do
    TUNA_WORK="$WORK/bench_m_work_m${M}"
    STDERR="$RESULTS/stderr_m${M}.txt"
    mkdir -p "$TUNA_WORK"

    "$TUNA" -k "$K" -m "$M" -t "$THREADS" -n "$N_PARTS" -dbg \
        -w "$TUNA_WORK/" "@$SUBSET_FOF" /dev/null \
        2>"$STDERR" || {
            echo "  m=$M FAILED — check $STDERR"
            continue
        }

    rm -rf "$TUNA_WORK"

    # Phase timings
    P1=$(grep  "^phase1:" "$STDERR" | awk -F: '{gsub(/s/,"",$2); printf "%.2f",$2}')
    P2=$(grep  "^phase2:" "$STDERR" | awk -F: '{gsub(/s/,"",$2); printf "%.2f",$2}')
    TOTAL=$(awk "BEGIN{printf \"%.2f\", $P1+$P2}")

    # Overflow stats
    OV_LINE=$(grep "k-mers went to overflow" "$STDERR" || echo "")
    if [ -n "$OV_LINE" ]; then
        OV_KMERS=$(echo "$OV_LINE" | awk '{print $2}')
        OV_PCT=$(echo  "$OV_LINE" | grep -oP '[\d.]+(?=% of total)')
    else
        OV_KMERS=0
        OV_PCT=0
    fi

    PARTS=$(grep "^tuna " "$STDERR" | grep -oP 'partitions=\K[0-9]+' || echo "$N_PARTS")

    # Resize stats (from -dbg output)
    N_RESIZES=$(grep "total resize events"  "$STDERR" | grep -oP ':\s*\K[0-9]+' || echo 0)
    RESIZE_S=$(grep  "total resize time"    "$STDERR" | grep -oP '[\d.]+(?=s)'   || echo 0)
    N_OV_RES=$(grep  "overflow-triggered"   "$STDERR" | grep -oP ':\s*\K[0-9]+' || echo 0)
    N_LD_RES=$(grep  "load-triggered"       "$STDERR" | grep -oP ':\s*\K[0-9]+' || echo 0)

    printf "%-4s  %-9s  %-9s  %-9s  %-10s  %-14s  %-9s  %-9s  %s/%s\n" \
        "$M" "$P1" "$P2" "$TOTAL" "$OV_PCT%" "$OV_KMERS" \
        "$N_RESIZES" "$RESIZE_S" "$N_OV_RES" "$N_LD_RES"

    echo "$M,$P1,$P2,$TOTAL,$OV_PCT,$OV_KMERS,$PARTS,$N_RESIZES,$RESIZE_S,$N_OV_RES,$N_LD_RES" >> "$CSV"
done

echo ""
echo "CSV: $CSV"
