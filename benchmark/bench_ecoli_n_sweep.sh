#!/usr/bin/env bash
# bench_ecoli_n_sweep.sh — tuna vs KMC across N E. coli genomes (log-scale N)
#
# Usage: bash bench_ecoli_n_sweep.sh [THREADS [MEM_GB]]
# Example: bash bench_ecoli_n_sweep.sh 4 12
#
# Path overrides via env vars:
#   TUNA=/path/to/tuna  KMC=/path/to/kmc  KMCDUMP=/path/to/kmc_dump
#   FOF=/path/to/ecoli_fof_3682.list  WORK=/path/to/workdir
#
# Output: $RESULTS/bench_n.csv
# Columns: n_files, tuna_p1_s, tuna_p2_s, tuna_total_s,
#          tuna_total_kmers, tuna_unique_kmers,
#          kmc_count_s, kmc_dump_s, kmc_total_s, kmc_unique_kmers

set -euo pipefail

THREADS=${1:-4}
MEM_GB=${2:-12}

# ── Paths ───────────────────────────────────────────────────────────────────
TUNA=${TUNA:-~/documents/giulio_colab/softs/tuna/build/tuna}
KMC=${KMC:-~/documents/giulio_colab/softs/kmc/bin/kmc}
KMCDUMP=${KMCDUMP:-~/documents/giulio_colab/softs/kmc/bin/kmc_dump}
FOF=${FOF:-~/tmp/data/ecoli_fof_3682.list}
WORK=${WORK:-~/tuna_bench_tmp}

# Expand ~ in paths
TUNA="${TUNA/#\~/$HOME}"
KMC="${KMC/#\~/$HOME}"
KMCDUMP="${KMCDUMP/#\~/$HOME}"
FOF="${FOF/#\~/$HOME}"
WORK="${WORK/#\~/$HOME}"

K=31
M=17   # best m for E. coli (lowest overflow, least resizes)

# Log-scale N values: 8 points from 1 to 3682
N_VALUES=(1 3 10 30 100 300 1000 3682)

RESULTS="$WORK/bench_n_sweep_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS"

CSV="$RESULTS/bench_n.csv"
echo "n_files,tuna_p1_s,tuna_p2_s,tuna_total_s,tuna_total_kmers,tuna_unique_kmers,kmc_count_s,kmc_dump_s,kmc_total_s,kmc_unique_kmers" > "$CSV"

echo "=== E. coli N-sweep: k=$K  m=$M  threads=$THREADS  mem=${MEM_GB}GB ==="
echo "Results: $RESULTS"
echo ""
printf "%-6s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-14s  %-14s\n" \
    "N" "tuna_p1" "tuna_p2" "tuna_tot" "kmc_cnt" "kmc_dmp" "kmc_tot" "tuna_total_K" "tuna_unique_K"
echo "──────  ─────────  ─────────  ─────────  ─────────  ─────────  ─────────  ──────────────  ──────────────"

for N in "${N_VALUES[@]}"; do

    # ── Subset fof ─────────────────────────────────────────────────────────
    SUBSET_FOF="$RESULTS/fof_${N}.list"
    head -n "$N" "$FOF" > "$SUBSET_FOF"

    # ─────────────────────────────────────────────────────────────────────────
    # tuna run
    # ─────────────────────────────────────────────────────────────────────────
    TUNA_WORK="$WORK/n_sweep_tuna_${N}"
    TUNA_ERR="$RESULTS/tuna_stderr_${N}.txt"
    mkdir -p "$TUNA_WORK"

    # Run without -hp so the "X k-mers processed, Y unique written" line appears.
    "$TUNA" -k "$K" -m "$M" -t "$THREADS" \
        -w "$TUNA_WORK/" "@$SUBSET_FOF" /dev/null \
        2>"$TUNA_ERR" || {
        echo "  N=$N tuna FAILED — see $TUNA_ERR"
        continue
    }
    rm -rf "$TUNA_WORK"

    P1=$(grep  "^phase1:" "$TUNA_ERR" | awk -F: '{gsub(/s/,"",$2); printf "%.3f",$2}')
    P2=$(grep  "^phase2:" "$TUNA_ERR" | awk -F: '{gsub(/s/,"",$2); printf "%.3f",$2}')
    TUNA_TOTAL=$(awk "BEGIN{printf \"%.3f\", $P1+$P2}")

    # Parse "N k-mers processed, M unique written"
    STATS_LINE=$(grep "k-mers processed" "$TUNA_ERR" || echo "")
    if [ -n "$STATS_LINE" ]; then
        TUNA_TOT_K=$(echo "$STATS_LINE" | awk '{print $1}')
        TUNA_UNI_K=$(echo "$STATS_LINE" | awk '{print $4}')
    else
        TUNA_TOT_K=0
        TUNA_UNI_K=0
    fi

    # ─────────────────────────────────────────────────────────────────────────
    # KMC run (count + dump)
    # ─────────────────────────────────────────────────────────────────────────
    KMC_WORK="$WORK/n_sweep_kmc_${N}"
    KMC_DB="$WORK/n_sweep_kmc_db_${N}"
    KMC_ERR="$RESULTS/kmc_stderr_${N}.txt"
    mkdir -p "$KMC_WORK"

    # Count
    T_KMC_CNT_START=$(date +%s%3N)
    "$KMC" -fm -k"$K" -t"$THREADS" -m"$MEM_GB" -ci1 \
        "@$SUBSET_FOF" "$KMC_DB" "$KMC_WORK" \
        >"$KMC_ERR" 2>&1 || {
        echo "  N=$N kmc count FAILED — see $KMC_ERR"
        rm -rf "$KMC_WORK" "${KMC_DB}.kmc_pre" "${KMC_DB}.kmc_suf" 2>/dev/null || true
        continue
    }
    T_KMC_CNT_END=$(date +%s%3N)
    KMC_CNT_S=$(awk "BEGIN{printf \"%.3f\", ($T_KMC_CNT_END - $T_KMC_CNT_START)/1000.0}")

    # Dump (to stdout, count unique k-mers with wc -l; suppress text output)
    T_KMC_DMP_START=$(date +%s%3N)
    KMC_UNI_K=$("$KMCDUMP" "$KMC_DB" /dev/stdout 2>/dev/null | wc -l)
    T_KMC_DMP_END=$(date +%s%3N)
    KMC_DMP_S=$(awk "BEGIN{printf \"%.3f\", ($T_KMC_DMP_END - $T_KMC_DMP_START)/1000.0}")
    KMC_TOTAL=$(awk "BEGIN{printf \"%.3f\", $KMC_CNT_S + $KMC_DMP_S}")

    rm -rf "$KMC_WORK" "${KMC_DB}.kmc_pre" "${KMC_DB}.kmc_suf"

    # ─────────────────────────────────────────────────────────────────────────
    # Print + record
    # ─────────────────────────────────────────────────────────────────────────
    TUNA_TOT_KM=$(awk "BEGIN{printf \"%.1f\", $TUNA_TOT_K/1000000.0}")
    TUNA_UNI_KM=$(awk "BEGIN{printf \"%.1f\", $TUNA_UNI_K/1000000.0}")

    printf "%-6s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-14s  %-14s\n" \
        "$N" "$P1" "$P2" "$TUNA_TOTAL" \
        "$KMC_CNT_S" "$KMC_DMP_S" "$KMC_TOTAL" \
        "${TUNA_TOT_KM}M" "${TUNA_UNI_KM}M"

    echo "$N,$P1,$P2,$TUNA_TOTAL,$TUNA_TOT_K,$TUNA_UNI_K,$KMC_CNT_S,$KMC_DMP_S,$KMC_TOTAL,$KMC_UNI_K" >> "$CSV"

done

echo ""
echo "CSV: $CSV"
