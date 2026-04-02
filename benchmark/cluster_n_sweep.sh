#!/usr/bin/env bash
# cluster_n_sweep.sh — sweep -n values on the cluster, collect phase timings.
#
# Usage:
#   bash cluster_n_sweep.sh <label> <input> <threads> <runs> <n values...>
#
# Examples:
#   bash cluster_n_sweep.sh ecoli @/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list 8 3 \
#       128 256 512 1024 2048 4096 8192 16384 32768 65536 131072
#
#   bash cluster_n_sweep.sh human \
#       /WORKS/vlevallois/data/dataset_genome_human/data/HG00438.maternal.f1_assembly_v2_genbank.fa.gz \
#       8 3 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288
#
# Output: TSV to stdout — redirect to a file:
#   bash cluster_n_sweep.sh ecoli ... > /WORKS/vlevallois/test_tuna/n_sweep_ecoli.tsv
#
# Notes:
#   - Input: plain path = single FASTA/FASTQ/gz file.
#             @ prefix  = file-of-files (tuna fof mode).
#   - RSS is peak resident set size in MB (from /usr/bin/time -v).
#   - If perf is available and perf_event_paranoid <= 1, cache counters can
#     be added later; this script requires only tuna + /usr/bin/time.

set -euo pipefail

LABEL=${1:?usage: cluster_n_sweep.sh <label> <input> <threads> <runs> <n...>}
INPUT=${2:?}
THREADS=${3:-8}
RUNS=${4:-3}
shift 4
N_VALUES=("$@")

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
WBASE=/WORKS/vlevallois/test_tuna/work_nsweep

mkdir -p "$WBASE"

if [ ! -x "$TUNA" ]; then
    echo "[ERROR] tuna binary not found or not executable: $TUNA" >&2
    exit 1
fi

printf "dataset\tn\trun\tphase1\tphase2\ttotal\trss_mb\n"

for N in "${N_VALUES[@]}"; do
    WD="${WBASE}_${LABEL}_${N}"
    mkdir -p "$WD"

    for i in $(seq 1 "$RUNS"); do
        TIMEFILE=$(mktemp /tmp/tuna_time_XXXXXX)

        /usr/bin/time -v -o "$TIMEFILE" \
            "$TUNA" -k 31 -m 21 -n "$N" -t "$THREADS" \
            "$INPUT" /dev/null -hp -w "$WD" \
            2>"$WD/run_${i}.stderr"

        P1=$(awk '/^phase1:/{gsub(/s/,"",$2); print $2}' "$WD/run_${i}.stderr")
        P2=$(awk '/^phase2:/{gsub(/s/,"",$2); print $2}' "$WD/run_${i}.stderr")
        TOTAL=$(awk "BEGIN{printf \"%.3f\", ${P1:-0} + ${P2:-0}}")
        RSS=$(grep "Maximum resident" "$TIMEFILE" | awk '{print $NF}')
        RSS_MB=$(awk "BEGIN{printf \"%.0f\", ${RSS:-0}/1024}")
        rm -f "$TIMEFILE"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$LABEL" "$N" "$i" "$P1" "$P2" "$TOTAL" "$RSS_MB"

        echo "[${LABEL}] n=${N} run=${i}: phase1=${P1}s phase2=${P2}s total=${TOTAL}s rss=${RSS_MB}MB" >&2
    done
done
