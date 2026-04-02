#!/usr/bin/env bash
# cluster_run_all.sh — run n-sweep for all datasets on the cluster.
#
# Usage: bash cluster_run_all.sh [threads] [runs]
#   threads — worker threads for tuna  (default: 8)
#   runs    — repetitions per (n, dataset) (default: 3)
#
# Output: /WORKS/vlevallois/test_tuna/n_sweep_<label>.tsv per dataset.
#
# Edit N_ECOLI / N_HUMAN below to change the partition values swept.
# Run datasets sequentially to avoid resource contention.

set -euo pipefail

THREADS=${1:-8}
RUNS=${2:-3}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$HERE/cluster_n_sweep.sh"
OUT=/WORKS/vlevallois/test_tuna

mkdir -p "$OUT"

# ── n values to sweep ─────────────────────────────────────────────────────────
# E. coli: tables are tiny; only interesting up to ~64K (TLB pressure appears).
N_ECOLI=(128 256 512 1024 2048 4096 8192 16384 32768 65536 131072)

# Human (single gz): full range including values beyond local sweet-spot.
# The cluster has 48MB L3/socket and more threads — optimal n will be higher.
N_HUMAN=(512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288)

# ── Dataset paths ─────────────────────────────────────────────────────────────
# Adjust if your cluster copies use different filenames.
ECOLI_FOF="/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list"
HUMAN_INPUT="/WORKS/vlevallois/data/dataset_genome_human/data/HG00438.maternal.f1_assembly_v2_genbank.fa.gz"

# ── E. coli sweep ─────────────────────────────────────────────────────────────
# Use first 200 files only (consistent with local benchmarks).
ECOLI_FOF_200="$OUT/ecoli_fof_200.list"
head -200 "$ECOLI_FOF" > "$ECOLI_FOF_200"
ECOLI_INPUT="@${ECOLI_FOF_200}"

echo "=== E. coli sweep  (t=$THREADS, ${#N_ECOLI[@]} n-values × $RUNS runs, 200 files) ===" >&2
bash "$SCRIPT" ecoli "$ECOLI_INPUT" "$THREADS" "$RUNS" "${N_ECOLI[@]}" \
    > "$OUT/n_sweep_ecoli.tsv"
echo "Saved: $OUT/n_sweep_ecoli.tsv" >&2

# ── Human sweep ───────────────────────────────────────────────────────────────
echo "=== Human sweep  (t=$THREADS, ${#N_HUMAN[@]} n-values × $RUNS runs) ===" >&2
bash "$SCRIPT" human "$HUMAN_INPUT" "$THREADS" "$RUNS" "${N_HUMAN[@]}" \
    > "$OUT/n_sweep_human.tsv"
echo "Saved: $OUT/n_sweep_human.tsv" >&2

echo "" >&2
echo "All done. Bring back results:" >&2
echo "  scp <cluster>:$OUT/n_sweep_ecoli.tsv ~/tuna_bench_tmp/" >&2
echo "  scp <cluster>:$OUT/n_sweep_human.tsv  ~/tuna_bench_tmp/" >&2
