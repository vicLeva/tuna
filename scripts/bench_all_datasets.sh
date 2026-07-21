#!/usr/bin/env bash
# scripts/bench_all_datasets.sh — tuna vs KMC, per-file sweep, 5 datasets
#
# Tools    : tuna · KMC
# Datasets : ecoli, human, salmonella, gut, tara (5 datasets)
# m        : 21 (fixed, all datasets)
# Mode     : count-only — output writing skipped entirely for both tools
#            (tuna: -co; KMC: -w), so wall-time reflects counting alone.
#
# Phase breakdown (both tools support it, verified against source):
#   tuna : own "phase1:"/"phase2:" stderr lines
#   KMC  : "1st stage:"/"2nd stage:" stderr lines — printed unconditionally
#          by kmc_CLI/kmc.cpp's print_summary() (undocumented in -h usage
#          text, but always emitted regardless of -w/-hp)
#
# Output: $RESULTS_DIR/bench_datasets.csv
#   dataset,file_idx,filename,tool,m,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_datasets/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"
: "${KMC:=""}"

THREADS=${THREADS:-8}
K=31
M=21
RAM_GB=256      # tuna -ram budget; also used as KMC -m

EXPES_ROOT="/WORKS/vlevallois/expes_tuna"
RESULTS_DIR="$EXPES_ROOT/bench_all_datasets_$(date +%Y%m%d_%H%M%S)"
WORKDIR="$EXPES_ROOT"

# Dataset registry: "name:fof_path:kmc_format:max_files"
DATASETS=(
    "ecoli:/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list:-fm:100"
    "human:/WORKS/vlevallois/data/dataset_genome_human/fof.list:-fm:10"
    "salmonella:/WORKS/vlevallois/data/dataset_pangenome_salmonella/fof.list:-fm:100"
    "gut:/WORKS/vlevallois/data/dataset_metagenome_gut/fof.list:-fm:100"
    "tara:/WORKS/vlevallois/data/dataset_metagenome_tara/fof.list:-fq:10"
)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for var in TUNA KMC; do
    [[ -z "${!var}" ]]              && { echo "[error] $var is not set"; err=1; }
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] not executable: ${!var}"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_datasets"
CSV="$RESULTS_DIR/bench_datasets.csv"
echo "dataset,file_idx,filename,tool,m,wall_s,rss_mb,phase1_s,phase2_s,unique_kmers" > "$CSV"

echo "[bench] tuna vs KMC — per-file dataset sweep (count-only, m=$M)"
echo "[bench] TUNA=$TUNA"
echo "[bench] KMC=$KMC"
echo "[bench] k=$K  m=$M  threads=$THREADS  ram=${RAM_GB}GB"
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
# Per-tool run functions
# Args: file work timefile stderrfile [kmc_fmt]
# =============================================================================

_run_tuna() {
    local file="$1" work="$2" tf="$3" se="$4"
    # -co: skip output writing after counting
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" -hp -co \
        -w "$work/" "$file" /dev/null \
        > /dev/null 2>"$se"
}

_run_kmc() {
    local file="$1" work="$2" tf="$3" se="$4" fmt="$5"
    # -w: without output — binary db is never written (true count-only)
    # Note: this cluster's KMC build prints its 1st/2nd-stage summary to
    # stdout, not stderr (differs from a newer git build tested locally) —
    # capture both streams into $se so se_val can find it.
    mkdir -p "$work/tmp"
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k${K} -ci1 -cs4294967295 "$fmt" -m${RAM_GB} -hp -t${THREADS} -w \
        "$file" "$work/out" "$work/tmp" \
        > "$se" 2>&1
    rm -rf "$work/tmp"
}

# =============================================================================
# Dispatcher: run one tool on one file and record results
# =============================================================================
run_one() {
    local tool="$1" ds="$2" file="$3" fidx="$4" fname="$5" kmc_fmt="$6"
    local tag="${ds}_f$(printf '%04d' "$fidx")_${tool}"
    local tf="$RESULTS_DIR/aux_datasets/${tag}.time"
    local se="$RESULTS_DIR/aux_datasets/${tag}.stderr"
    local work="$WORKDIR/ds_work_${tag}"
    mkdir -p "$work"

    local ok=true
    case "$tool" in
        tuna) _run_tuna "$file" "$work" "$tf" "$se" ;;
        kmc)  _run_kmc  "$file" "$work" "$tf" "$se" "$kmc_fmt" ;;
    esac || ok=false

    rm -rf "$work"
    $ok || { echo "    [FAIL] $tool"; return; }

    local wall rss p1 p2 unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    case "$tool" in
        tuna) p1=$(se_val "phase1"    "$se")
              p2=$(se_val "phase2"    "$se")
              unique=$(se_val "unique_kmers" "$se") ;;
        kmc)  p1=$(se_val "1st stage" "$se")
              p2=$(se_val "2nd stage" "$se")
              unique=$(grep "No. of unique k-mers" "$se" | awk -F: '{gsub(/^ +/,"",$2); print $2}') ;;
    esac

    printf "    [%-4s] wall=%8ss  p1=%8ss  p2=%8ss  RSS=%6sMB  unique=%s\n" \
        "$tool" "$wall" "$p1" "$p2" "$rss" "$unique"
    echo "$ds,$fidx,$fname,$tool,$M,$wall,$rss,$p1,$p2,$unique" >> "$CSV"
}

# =============================================================================
# Main loop
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"

for ENTRY in "${DATASETS[@]}"; do
    IFS=: read -r DS_NAME FOF KMC_FMT DS_MAX <<< "$ENTRY"

    if [[ ! -f "$FOF" ]]; then
        echo "  [SKIP] $DS_NAME: fof not found: $FOF"
        continue
    fi

    TOTAL=$(wc -l < "$FOF")
    N=$(( TOTAL < DS_MAX ? TOTAL : DS_MAX ))
    mapfile -t FILES < <(head -n "$N" "$FOF")

    echo ""
    echo "──── $DS_NAME  ($N / $TOTAL files) ────────────────────"

    for IDX in "${!FILES[@]}"; do
        FILE="${FILES[$IDX]}"
        FNAME=$(basename "$FILE")
        FIDX=$(( IDX + 1 ))
        printf "  [%d/%d] %s\n" "$FIDX" "$N" "$FNAME"
        # Alternate order each file so neither tool consistently gets warm cache.
        if (( FIDX % 2 == 1 )); then
            run_one kmc  "$DS_NAME" "$FILE" "$FIDX" "$FNAME" "$KMC_FMT"
            run_one tuna "$DS_NAME" "$FILE" "$FIDX" "$FNAME" ""
        else
            run_one tuna "$DS_NAME" "$FILE" "$FIDX" "$FNAME" ""
            run_one kmc  "$DS_NAME" "$FILE" "$FIDX" "$FNAME" "$KMC_FMT"
        fi
    done
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
