#!/usr/bin/env bash
# scripts/bench_all_datasets_dump.sh — tuna vs KMC, per-file sweep, 5 datasets
#
# Bis of bench_all_datasets.sh with real output written to disk instead of
# count-only, to measure the actual dump cost (vs. new_runs_v2's count-only
# baseline).
#
# Tools    : tuna · KMC
# Datasets : ecoli, human, salmonella, gut, tara (5 datasets)
# m        : 21 (fixed, all datasets)
# Mode     : full dump — tuna writes its TSV, KMC writes its binary db then
#            dumps it to ASCII via kmc_dump. Output files are deleted right
#            after each run's timing is captured (not kept — full dumps can
#            be very large, e.g. ~100GB for a human genome TSV).
#
# Phase breakdown:
#   tuna : own "phase1:"/"phase2:" stderr lines (now includes output-write
#          cost interleaved into phase2 — tuna writes per-partition, so the
#          dump cost can't be isolated directly; compare phase2_s here
#          against new_runs_v2's count-only phase2_s to estimate it)
#   KMC  : "1st stage:"/"2nd stage:" stderr lines, plus a separate timed
#          kmc_dump step (dump_s) — KMC's own binary-write cost is folded
#          into its stage times same as tuna; compare against new_runs_v2
#          to estimate it. dump_s here is the extra ASCII-dump step only.
#
# Output: $RESULTS_DIR/bench_datasets_dump.csv
#   dataset,file_idx,filename,tool,m,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_datasets/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"
: "${KMC:=""}"
: "${KMC_DUMP:=""}"

THREADS=${THREADS:-8}
K=31
M=21
RAM_GB=256      # tuna -ram budget; also used as KMC -m

EXPES_ROOT="/WORKS/vlevallois/expes_tuna"
RESULTS_DIR="$EXPES_ROOT/bench_all_datasets_dump_$(date +%Y%m%d_%H%M%S)"
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
for var in TUNA KMC KMC_DUMP; do
    [[ -z "${!var}" ]]              && { echo "[error] $var is not set"; err=1; }
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] not executable: ${!var}"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_datasets"
CSV="$RESULTS_DIR/bench_datasets_dump.csv"
echo "dataset,file_idx,filename,tool,m,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers" > "$CSV"

echo "[bench] tuna vs KMC — per-file dataset sweep (full dump, m=$M)"
echo "[bench] TUNA=$TUNA"
echo "[bench] KMC=$KMC"
echo "[bench] KMC_DUMP=$KMC_DUMP"
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
    # Full dump: real TSV output instead of /dev/null, no -co.
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" -hp \
        -w "$work/" "$file" "$work/out.tsv" \
        > /dev/null 2>"$se"
    local rc=$?
    rm -f "$work/out.tsv"
    return "$rc"
}

_run_kmc() {
    local file="$1" work="$2" tf="$3" se="$4" fmt="$5" dt="$6"
    # Full dump: real binary db (no -w), then a separately-timed kmc_dump
    # step to get ASCII output. Note: this cluster's KMC build prints its
    # 1st/2nd-stage summary to stdout, not stderr — capture both streams.
    mkdir -p "$work/tmp"
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k${K} -ci1 -cs4294967295 "$fmt" -m${RAM_GB} -hp -t${THREADS} \
        "$file" "$work/out" "$work/tmp" \
        > "$se" 2>&1
    local rc=$?
    rm -rf "$work/tmp"
    if [[ "$rc" -ne 0 ]]; then return "$rc"; fi

    /usr/bin/time -v -o "$dt" \
        "$KMC_DUMP" "$work/out" "$work/out.tsv" \
        > /dev/null 2>>"$se"
    rc=$?
    rm -f "$work/out.tsv" "$work/out.kmc_pre" "$work/out.kmc_suf"
    return "$rc"
}

# =============================================================================
# Dispatcher: run one tool on one file and record results
# =============================================================================
run_one() {
    local tool="$1" ds="$2" file="$3" fidx="$4" fname="$5" kmc_fmt="$6"
    local tag="${ds}_f$(printf '%04d' "$fidx")_${tool}"
    local tf="$RESULTS_DIR/aux_datasets/${tag}.time"
    local se="$RESULTS_DIR/aux_datasets/${tag}.stderr"
    local dt="$RESULTS_DIR/aux_datasets/${tag}.dumptime"
    local work="$WORKDIR/ds_work_${tag}"
    mkdir -p "$work"

    local rc=0
    case "$tool" in
        tuna) _run_tuna "$file" "$work" "$tf" "$se"; rc=$? ;;
        kmc)  _run_kmc  "$file" "$work" "$tf" "$se" "$kmc_fmt" "$dt"; rc=$? ;;
    esac

    rm -rf "$work"
    [[ "$rc" -eq 0 ]] || { echo "    [FAIL] $tool"; return; }

    local wall rss p1 p2 dump unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    case "$tool" in
        tuna) p1=$(se_val "phase1"    "$se")
              p2=$(se_val "phase2"    "$se")
              dump=na
              unique=$(se_val "unique_kmers" "$se") ;;
        kmc)  p1=$(se_val "1st stage" "$se")
              p2=$(se_val "2nd stage" "$se")
              dump=$(wall_from_file "$dt")
              unique=$(grep "No. of unique k-mers" "$se" | awk -F: '{gsub(/^ +/,"",$2); print $2}') ;;
    esac

    printf "    [%-4s] wall=%8ss  p1=%8ss  p2=%8ss  dump=%8ss  RSS=%6sMB  unique=%s\n" \
        "$tool" "$wall" "$p1" "$p2" "$dump" "$rss" "$unique"
    echo "$ds,$fidx,$fname,$tool,$M,$wall,$rss,$p1,$p2,$dump,$unique" >> "$CSV"
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
