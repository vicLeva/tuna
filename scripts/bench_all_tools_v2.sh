#!/usr/bin/env bash
# scripts/bench_all_tools_v2.sh — Rob's tuna, assembly datasets, thread sweep
#
# Mirrors bench_all_tools.sh but tuna-only (Rob's fork).
# Datasets : first genome from ECOLI_FOF and HUMAN_FOF
# Threads  : configurable sweep (default: 1 2 4 6 8 10 16 22 28 32)
# m        : 21 for E. coli, 23 for Human (assembled genome, highly repetitive)
# n        : auto-tuned (Rob's L3-aware logic; not overridden)
#
# Output: $RESULTS_DIR/bench_tools.csv
#   dataset,file,threads,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_tools/

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

K=31
RAM_GB=256
THREADS_LIST=(1 2 4 6 8 10 16 22 28 32)

# m per dataset: human assemblies need m=23 (highly repetitive, large genome)
declare -A DS_M=([ecoli]=21 [human]=23)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
[[ -z "$TUNA" ]] && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
[[ ! -f "$ECOLI_FOF" ]] && { echo "[error] not found: $ECOLI_FOF"; err=1; }
[[ ! -f "$HUMAN_FOF" ]] && { echo "[error] not found: $HUMAN_FOF"; err=1; }
[[ "$err" -eq 1 ]] && exit 1

ECOLI_GENOME=$(head -1 "$ECOLI_FOF")
HUMAN_GENOME=$(head -1 "$HUMAN_FOF")
[[ ! -f "$ECOLI_GENOME" ]] && { echo "[error] E. coli genome not found: $ECOLI_GENOME"; exit 1; }
[[ ! -f "$HUMAN_GENOME" ]] && { echo "[error] Human genome not found: $HUMAN_GENOME"; exit 1; }

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_tools"
CSV="$RESULTS_DIR/bench_tools.csv"
echo "dataset,file,threads,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers" > "$CSV"

echo "[bench] Rob's tuna — assembly thread sweep"
echo "[bench] TUNA=$TUNA"
echo "[bench] k=$K  ram=${RAM_GB}GB  threads=${THREADS_LIST[*]}"
echo "[bench] E. coli: $(basename "$ECOLI_GENOME")"
echo "[bench] Human  : $(basename "$HUMAN_GENOME")"
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
# Extract a "KEY: VALUEs" line from tuna stderr (strips trailing 's' if present)
se_val() { grep "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2; gsub(/^ +| +s$/,"",v); print v}' || echo na; }

# =============================================================================
# Run function
# =============================================================================
run_tuna() {
    local ds="$1" genome="$2" threads="$3"
    local m="${DS_M[$ds]}"
    local tag="${ds}_t${threads}"
    local tf="$RESULTS_DIR/aux_tools/${tag}.time"
    local se="$RESULTS_DIR/aux_tools/${tag}.stderr"
    local out="$WORKDIR/tuna_tools_${tag}.tsv"
    local work="$WORKDIR/tuna_tools_work_${tag}"
    mkdir -p "$work"

    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$m" -t "$threads" -ram "$RAM_GB" -hp \
        -w "$work/" "$genome" "$out" \
        > /dev/null 2>"$se" \
    || { echo "  [FAIL] tuna $ds t=$threads"; rm -rf "$work" "$out"; return; }
    rm -rf "$work" "$out"

    local wall rss p1 p2 n_parts unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")
    n_parts=$(se_val "n_parts" "$se")
    unique=$(se_val "unique_kmers" "$se")

    printf "  %-6s t=%-3d m=%d  wall=%8.3fs  p1=%ss p2=%ss  RSS=%s MB  n=%s\n" \
        "$ds" "$threads" "$m" "$wall" "$p1" "$p2" "$rss" "$n_parts"
    echo "$ds,$(basename "$genome"),$threads,$m,$wall,$rss,$p1,$p2,$n_parts,$unique" >> "$CSV"
}

# =============================================================================
# Main loop
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"
for threads in "${THREADS_LIST[@]}"; do
    echo "── t=$threads ──────────────────────────────────────────"
    run_tuna ecoli "$ECOLI_GENOME" "$threads"
    run_tuna human "$HUMAN_GENOME" "$threads"
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
