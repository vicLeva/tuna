#!/usr/bin/env bash
# perf_cache_bench.sh — LLC cache-miss analysis vs n_parts
#
# Usage: bash perf_cache_bench.sh <input_file> [THREADS] [K] [M] [PHASE]
#
#   input_file — single FASTA/FASTQ (plain or .gz)
#   THREADS    — default 4
#   K          — default 31
#   M          — default 21
#   PHASE      — "1" (partition only, -tp), "2" (full), "both" (default)
#
# Sweeps N_VALUES powers-of-2 from 512 to 32768.
# For each n_parts:
#   • [phase=both] runs full pipeline, reports phase1+phase2 LLC misses
#   • [phase=1]    runs with -tp (partition only), phase1 LLC misses
#   • [phase=2]    runs full, then subtracts -tp run for phase2 estimate
#
# Output: $OUT_DIR/cache_bench.csv
#   n_parts, phase, llc_load_misses, llc_store_misses, superkmers, kmers,
#   misses_per_superkmer, misses_per_kmer, phase_wall_s
#
# Requirements: perf, /usr/bin/time -v, tuna binary
#
# Expected result: after n_parts ≤ 8192 cap, LLC-load-misses/superkmer
# should be roughly constant (≈1–2 per superkmer for phase2 hash table
# probes), with a sharp rise above 8192 (writer-buffer cache thrashing).

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────

INPUT=${1:?Usage: $0 <input_file> [THREADS] [K] [M] [PHASE]}
THREADS=${2:-4}
K=${3:-31}
M=${4:-21}
PHASE=${5:-both}   # "1", "2", or "both"

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
WORK=/WORKS/vlevallois/test_tuna
OUT_DIR="$WORK/cache_bench_$(date +%Y%m%d_%H%M%S)"

N_VALUES=(512 1024 2048 4096 8192 16384 32768)

mkdir -p "$OUT_DIR"

# ── Helpers ───────────────────────────────────────────────────────────────────

# Extract a perf stat event value (removes commas, returns integer)
# perf stat writes to stderr; we capture it via file descriptor trick.
parse_perf() {
    local event="$1" perf_file="$2"
    grep -E "^[[:space:]]*[0-9,].*${event}" "$perf_file" \
        | head -1 | awk '{gsub(/,/,""); print $1}'
}

parse_tuna_field() {
    local field="$1" stderr_file="$2"
    # field like "phase1", "phase2", "superkmers"
    grep "^${field}:" "$stderr_file" | awk -F: '{gsub(/s/,"",$2); printf "%.3f", $2}'
}

run_once() {
    local n="$1" extra_flags="$2" tag="$3"
    local tuna_work="$WORK/tuna_work_cache_n${n}_${tag}"
    local tuna_out="$WORK/tuna_out_cache_$$_${tag}.tsv"
    local tuna_stderr="$OUT_DIR/n${n}_${tag}.stderr"
    local perf_out="$OUT_DIR/n${n}_${tag}.perf"

    mkdir -p "$tuna_work"

    # Use -o to write perf output directly to a file, independent of stderr.
    perf stat -o "$perf_out" \
        -e LLC-load-misses,LLC-store-misses,cache-misses \
        -- \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -n "$n" -hp \
            $extra_flags \
            -w "$tuna_work/" "$INPUT" "$tuna_out" \
        2>"$tuna_stderr"

    rm -rf "$tuna_work" "$tuna_out" 2>/dev/null || true
}

# ── CSV header ────────────────────────────────────────────────────────────────

CSV="$OUT_DIR/cache_bench.csv"
echo "n_parts,phase,llc_load_misses,llc_store_misses,superkmers,kmers,misses_per_sk,misses_per_kmer,wall_s" \
    > "$CSV"

echo "=== LLC cache-miss sweep ==="
echo "    input:   $INPUT"
echo "    k=$K  m=$M  threads=$THREADS"
echo "    phase:   $PHASE"
echo "    n sweep: ${N_VALUES[*]}"
echo "    output:  $OUT_DIR"
echo ""
printf "%-8s  %-6s  %-14s  %-14s  %-12s  %-12s  %-14s  %-12s\n" \
    "n_parts" "phase" "LLC-load-miss" "LLC-stor-miss" "superkmers" "kmers" "miss/superkmer" "wall_s"
printf '%s\n' "─────────────────────────────────────────────────────────────────────────────────────────────────"

# ── Main sweep ────────────────────────────────────────────────────────────────

for N in "${N_VALUES[@]}"; do

    if [ "$PHASE" = "1" ] || [ "$PHASE" = "both" ]; then

        run_once "$N" "-tp" "p1"
        p1_llc_load=$(parse_perf "LLC-load-misses"  "$OUT_DIR/n${N}_p1.perf")
        p1_llc_stor=$(parse_perf "LLC-store-misses" "$OUT_DIR/n${N}_p1.perf")
        p1_sk=$(parse_tuna_field "superkmers" "$OUT_DIR/n${N}_p1.stderr")
        p1_km=$(grep "^phase1:" "$OUT_DIR/n${N}_p1.stderr" | wc -l)  # placeholder — kmers not in stderr for -tp
        # Actually for -tp, kmers are logged via progress line — extract from non-structured output
        # Use superkmers as proxy; k-mers ≈ sk × avg_sk_len
        p1_wall=$(parse_tuna_field "phase1" "$OUT_DIR/n${N}_p1.stderr")

        # kmers: parse from the progress line "N sequences, M k-mers"
        p1_km=$(grep "k-mers" "$OUT_DIR/n${N}_p1.stderr" | awk '{gsub(/,/,""); for(i=1;i<=NF;i++) if($i~/^[0-9]+$/ && $(i+1)=="k-mers") print $i}' | head -1)
        p1_km=${p1_km:-0}
        p1_sk=${p1_sk:-0}
        p1_llc_load=${p1_llc_load:-0}
        p1_llc_stor=${p1_llc_stor:-0}
        p1_wall=${p1_wall:-0}

        if [ "$p1_sk" -gt 0 ] 2>/dev/null; then
            p1_miss_sk=$(awk "BEGIN{printf \"%.3f\", ($p1_llc_load+$p1_llc_stor)/$p1_sk}")
            p1_miss_km=$(awk "BEGIN{printf \"%.4f\", ($p1_llc_load+$p1_llc_stor)/($p1_km>0?$p1_km:1)}")
        else
            p1_miss_sk="N/A"; p1_miss_km="N/A"
        fi

        printf "%-8s  %-6s  %-14s  %-14s  %-12s  %-12s  %-14s  %-12s\n" \
            "$N" "phase1" "$p1_llc_load" "$p1_llc_stor" "$p1_sk" "$p1_km" "$p1_miss_sk" "$p1_wall"
        echo "${N},phase1,${p1_llc_load},${p1_llc_stor},${p1_sk},${p1_km},${p1_miss_sk},${p1_miss_km},${p1_wall}" >> "$CSV"
    fi

    if [ "$PHASE" = "2" ] || [ "$PHASE" = "both" ]; then

        run_once "$N" "" "full"
        f_llc_load=$(parse_perf "LLC-load-misses"  "$OUT_DIR/n${N}_full.perf")
        f_llc_stor=$(parse_perf "LLC-store-misses" "$OUT_DIR/n${N}_full.perf")
        f_sk=$(parse_tuna_field "superkmers" "$OUT_DIR/n${N}_full.stderr")
        f_wall_p2=$(parse_tuna_field "phase2"  "$OUT_DIR/n${N}_full.stderr")
        f_km=$(grep "k-mers" "$OUT_DIR/n${N}_full.stderr" | awk '{gsub(/,/,""); for(i=1;i<=NF;i++) if($i~/^[0-9]+$/ && $(i+1)=="k-mers") print $i}' | head -1)
        f_km=${f_km:-0}
        f_sk=${f_sk:-0}

        # If we also ran phase1, subtract its LLC misses to isolate phase2
        if [ "$PHASE" = "both" ] && [ -f "$OUT_DIR/n${N}_p1.perf" ]; then
            p2_llc_load=$(awk -v a="${f_llc_load:-0}" -v b="${p1_llc_load:-0}" 'BEGIN{x=a-b; print (x<0?0:x)}')
            p2_llc_stor=$(awk -v a="${f_llc_stor:-0}" -v b="${p1_llc_stor:-0}" 'BEGIN{x=a-b; print (x<0?0:x)}')
        else
            p2_llc_load=${f_llc_load:-0}
            p2_llc_stor=${f_llc_stor:-0}
        fi

        if [ "${f_sk:-0}" -gt 0 ] 2>/dev/null; then
            p2_miss_sk=$(awk "BEGIN{printf \"%.3f\", ($p2_llc_load+$p2_llc_stor)/($f_sk>0?$f_sk:1)}")
            p2_miss_km=$(awk "BEGIN{printf \"%.4f\", ($p2_llc_load+$p2_llc_stor)/($f_km>0?$f_km:1)}")
        else
            p2_miss_sk="N/A"; p2_miss_km="N/A"
        fi

        printf "%-8s  %-6s  %-14s  %-14s  %-12s  %-12s  %-14s  %-12s\n" \
            "$N" "phase2" "$p2_llc_load" "$p2_llc_stor" "$f_sk" "$f_km" "$p2_miss_sk" "$f_wall_p2"
        echo "${N},phase2,${p2_llc_load},${p2_llc_stor},${f_sk},${f_km},${p2_miss_sk},${p2_miss_km},${f_wall_p2}" >> "$CSV"
    fi

done

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "=== Done ==="
echo "Results: $CSV  ($(wc -l < "$CSV") rows)"
echo ""
echo "Expected: phase1 miss/superkmer flat ≤8192, rising above (writer-buffer thrashing)"
echo "          phase2 miss/superkmer ≈1–2 (one LLC miss per hash table probe, hidden by prefetch)"
echo ""
column -t -s, "$CSV"
