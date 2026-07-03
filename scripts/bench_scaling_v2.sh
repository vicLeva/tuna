#!/usr/bin/env bash
# scripts/bench_scaling_v2.sh — Rob's tuna, n-files scaling sweep
#
# Mirrors bench_scaling.sh but tuna-only (Rob's fork).
# Same N values as original for direct comparison with 3-tool results.
#
# Datasets:
#   ecoli  : N = 1 2 3 5 10 20 50 100 200 500 1000 1500 2000 2500 3000 3500
#   human  : N = 1 2 3 4 5 6 7 8 9 10 15 20 25 30
#   gallus : N = 1 2 3 5 8 10 15 (all files)     [reads, new]
#   human3 : N = 1 2 3 5 10 15 20 25 30 36 (all) [reads, new]
#
# m : 21 for ecoli/gallus/human3, 23 for human assemblies
# n : auto-tuned (Rob's L3-aware logic)
#
# Output: $RESULTS_DIR/bench_scaling.csv
#   dataset,n_files,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers
# Aux files in $RESULTS_DIR/aux_scaling/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"

RESULTS_DIR="/WORKS/vlevallois/expes_tuna/rob_v2"
WORKDIR="/WORKS/vlevallois/expes_tuna"

ECOLI_FOF="/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list"
HUMAN_FOF="/WORKS/vlevallois/data/dataset_genome_human/fof.list"
GALLUS_FOF="/WORKS/vlevallois/data/dataset_reads_gallus/fof.list"
HUMAN3_FOF="/WORKS/vlevallois/data/dataset_reads_human3/fof.list"

K=31
RAM_GB=256
THREADS=${THREADS:-8}

ECOLI_NS=(1 2 3 5 10 20 50 100 200 500 1000 1500 2000 2500 3000 3500)
HUMAN_NS=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)
GALLUS_NS=(1 2 3 5 8 10 15)
HUMAN3_NS=(1 2 3 5 10 15 20 25 30 36)

# m per dataset
declare -A DS_M=([ecoli]=21 [human]=23 [gallus]=21 [human3]=21)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
[[ -z "$TUNA" ]] && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
for var in ECOLI_FOF HUMAN_FOF; do
    [[ ! -f "${!var}" ]] && { echo "[error] not found: ${!var}"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

for var in GALLUS_FOF HUMAN3_FOF; do
    [[ ! -f "${!var}" ]] && echo "[warn] not found (skipping): ${!var}"
done

ECOLI_TOTAL=$(wc -l < "$ECOLI_FOF")
HUMAN_TOTAL=$(wc -l < "$HUMAN_FOF")
GALLUS_TOTAL=$([[ -f "$GALLUS_FOF" ]] && wc -l < "$GALLUS_FOF" || echo 0)
HUMAN3_TOTAL=$([[ -f "$HUMAN3_FOF" ]] && wc -l < "$HUMAN3_FOF" || echo 0)

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_scaling"
CSV="$RESULTS_DIR/bench_scaling.csv"
echo "dataset,n_files,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers" > "$CSV"

echo "[bench] Rob's tuna — n-files scaling sweep"
echo "[bench] TUNA=$TUNA"
echo "[bench] k=$K  threads=$THREADS  ram=${RAM_GB}GB"
echo "[bench] E. coli fof:  $ECOLI_TOTAL files"
echo "[bench] Human fof:    $HUMAN_TOTAL files"
echo "[bench] Gallus fof:   $GALLUS_TOTAL files"
echo "[bench] Human3 fof:   $HUMAN3_TOTAL files"
echo "[bench] Results: $CSV"

# =============================================================================
# Helpers
# =============================================================================
wall_from_file() {
    local t; t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    awk -F: '{if(NF==3) printf "%.3f",$1*3600+$2*60+$3; else printf "%.3f",$1*60+$2}' <<< "$t"
}
rss_mb() {
    local kb; kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}
se_val() { grep "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2; gsub(/^ +| +s$/,"",v); print v}' || echo na; }

# =============================================================================
# Run function
# =============================================================================
run_tuna() {
    local ds="$1" n="$2" fof="$3"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then
        echo "  [skip] $ds n=$n > fof ($total)"
        return
    fi

    local m="${DS_M[$ds]}"
    local subfof="$RESULTS_DIR/aux_scaling/subfof_${ds}_n${n}.list"
    head -n "$n" "$fof" > "$subfof"

    local tag="${ds}_n${n}"
    local tf="$RESULTS_DIR/aux_scaling/${tag}.time"
    local se="$RESULTS_DIR/aux_scaling/${tag}.stderr"
    local out="$WORKDIR/tuna_scaling_${tag}.tsv"
    local work="$WORKDIR/tuna_scaling_work_${tag}"
    mkdir -p "$work"

    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$m" -t "$THREADS" -ram "$RAM_GB" -hp \
        -w "$work/" "@$subfof" "$out" \
        > /dev/null 2>"$se" \
    || { echo "  [FAIL] tuna $ds n=$n"; rm -rf "$work" "$subfof" "$out"; return; }
    rm -rf "$work" "$subfof" "$out"

    local wall rss p1 p2 n_parts unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")
    n_parts=$(se_val "n_parts" "$se")
    unique=$(se_val "unique_kmers" "$se")

    printf "  %-8s n=%-5d m=%d  wall=%8.3fs  p1=%ss p2=%ss  RSS=%s MB  n=%s\n" \
        "$ds" "$n" "$m" "$wall" "$p1" "$p2" "$rss" "$n_parts"
    echo "$ds,$n,$wall,$rss,$p1,$p2,$n_parts,$unique" >> "$CSV"
}

# =============================================================================
# Main loop: full dataset before moving to next (mirrors bench_scaling.sh order)
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"

echo "════ ecoli ══════════════════════════════════════════════"
for n in "${ECOLI_NS[@]}"; do
    run_tuna ecoli "$n" "$ECOLI_FOF"
done

echo "════ human ══════════════════════════════════════════════"
for n in "${HUMAN_NS[@]}"; do
    run_tuna human "$n" "$HUMAN_FOF"
done

if [[ -f "$GALLUS_FOF" && "$GALLUS_TOTAL" -gt 0 ]]; then
    echo "════ gallus (reads) ══════════════════════════════════"
    for n in "${GALLUS_NS[@]}"; do
        run_tuna gallus "$n" "$GALLUS_FOF"
    done
fi

if [[ -f "$HUMAN3_FOF" && "$HUMAN3_TOTAL" -gt 0 ]]; then
    echo "════ human3 (reads) ══════════════════════════════════"
    for n in "${HUMAN3_NS[@]}"; do
        run_tuna human3 "$n" "$HUMAN3_FOF"
    done
fi

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
