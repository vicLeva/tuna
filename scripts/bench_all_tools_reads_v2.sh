#!/usr/bin/env bash
# scripts/bench_all_tools_reads_v2.sh — reads datasets, t=8, 2 RAM experiments
#
# Tools    : tuna_rob · KMC · FastK   (tuna dev excluded)
# Datasets : G. gallus reads (SRX043656) and H. sapiens reads (ERR174324-ERR174341)
#            Both passed as @fof (all files counted together).
# Threads  : fixed at 8 (no sweep)
# m        : 21 for both (short reads; m=23 is for assembled genomes)
# n        : auto-tuned for tuna_rob
#
# Experiments (RAM budget × tool subset):
#   1. RAM=5GB    : KMC, tuna_rob            (FastK excluded)
#   2. RAM=512GB  : KMC, tuna_rob, FastK
#
# Output skipped for all tools:
#   tuna_rob : output path is /dev/null; also passes -co
#   KMC      : binary db written then deleted; dump step omitted
#   FastK    : .ktab written then deleted; Tabex step omitted
#
# Output: $RESULTS_DIR/bench_reads.csv
#   dataset,n_files,threads,ram_gb,tool,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers
# Aux files (timefile, stderr) in $RESULTS_DIR/aux_reads/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA_ROB:=""}"  # tuna_rob binary
: "${KMC:=""}"       # kmc binary
: "${FASTK:=""}"     # FastK binary

RESULTS_DIR="/WORKS/vlevallois/expes_tuna/rob_v2"
WORKDIR="/WORKS/vlevallois/expes_tuna"

GALLUS_FOF="/WORKS/vlevallois/data/dataset_reads_gallus/fof.list"
HUMAN3_FOF="/WORKS/vlevallois/data/dataset_reads_human3/fof.list"

K=31
M=21
THREADS=8

# "ram_gb:tool1 tool2 ..."
EXPERIMENTS=(
    "5:kmc tuna_rob"
    "512:kmc tuna_rob fastk"
)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for var in TUNA_ROB KMC FASTK; do
    [[ -z "${!var}" ]]              && { echo "[error] $var is not set"; err=1; }
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] not executable: ${!var}"; err=1; }
done
[[ ! -f "$GALLUS_FOF" ]] && { echo "[error] not found: $GALLUS_FOF"; err=1; }
[[ ! -f "$HUMAN3_FOF" ]] && { echo "[error] not found: $HUMAN3_FOF"; err=1; }
[[ "$err" -eq 1 ]] && exit 1

GALLUS_N=$(wc -l < "$GALLUS_FOF")
HUMAN3_N=$(wc -l < "$HUMAN3_FOF")

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_reads"
CSV="$RESULTS_DIR/bench_reads.csv"
echo "dataset,n_files,threads,ram_gb,tool,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers" > "$CSV"

echo "[bench] reads benchmark — t=$THREADS, 2 RAM experiments"
echo "[bench] TUNA_ROB=$TUNA_ROB"
echo "[bench] KMC=$KMC"
echo "[bench] FASTK=$FASTK"
echo "[bench] k=$K  m=$M  threads=$THREADS"
for exp in "${EXPERIMENTS[@]}"; do
    echo "[bench] experiment: ram=${exp%%:*}GB  tools=${exp#*:}"
done
echo "[bench] G. gallus  fof: $GALLUS_FOF ($GALLUS_N files)"
echo "[bench] H. sapiens fof: $HUMAN3_FOF ($HUMAN3_N files)"
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
# Args: fof work timefile stderrfile threads ram_gb
# =============================================================================

_run_tuna_rob() {
    local fof="$1" work="$2" tf="$3" se="$4" threads="$5" ram_gb="$6"
    # -co: skip output writing after counting
    /usr/bin/time -v -o "$tf" \
        "$TUNA_ROB" -k "$K" -m "$M" -t "$threads" -ram "$ram_gb" -hp -co \
        -w "$work/" "@$fof" /dev/null \
        > /dev/null 2>"$se"
}

_run_kmc() {
    local fof="$1" work="$2" tf="$3" se="$4" threads="$5" ram_gb="$6"
    # -fq: FASTQ input; @fof: list of files; no dump step (binary db deleted after)
    mkdir -p "$work/tmp"
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k${K} -ci1 -cs4294967295 -fq -m${ram_gb} -hp -t${threads} \
        "@$fof" "$work/out" "$work/tmp" \
        > /dev/null 2>"$se"
    rm -f "$work/out"* && rm -rf "$work/tmp"
}

_run_fastk() {
    local fof="$1" work="$2" tf="$3" se="$4" threads="$5" ram_gb="$6"
    # .fastq.gz is recognised natively by FastK; no symlinks needed
    # No Tabex dump step; .ktab files deleted after
    mkdir -p "$work/tmp"
    mapfile -t files < "$fof"
    /usr/bin/time -v -o "$tf" \
        "$FASTK" -k${K} -t1 -T${threads} -M${ram_gb} \
        -N"$work/out" -P"$work/tmp" \
        "${files[@]}" \
        > /dev/null 2>"$se"
    rm -f "$work/out"* && rm -rf "$work/tmp"
}

# =============================================================================
# Dispatcher: run one tool on one dataset and record results
# =============================================================================
run_one() {
    local tool="$1" ds="$2" fof="$3" n_files="$4" threads="$5" ram_gb="$6"
    local tag="${ds}_${tool}_t${threads}_ram${ram_gb}"
    local tf="$RESULTS_DIR/aux_reads/${tag}.time"
    local se="$RESULTS_DIR/aux_reads/${tag}.stderr"
    local work="$WORKDIR/reads_work_${tag}"
    mkdir -p "$work"

    local ok=true
    case "$tool" in
        tuna_rob) _run_tuna_rob "$fof" "$work" "$tf" "$se" "$threads" "$ram_gb" ;;
        kmc)      _run_kmc      "$fof" "$work" "$tf" "$se" "$threads" "$ram_gb" ;;
        fastk)    _run_fastk    "$fof" "$work" "$tf" "$se" "$threads" "$ram_gb" ;;
    esac || ok=false

    rm -rf "$work"
    $ok || { echo "  [FAIL] $tool $ds t=$threads ram=${ram_gb}GB"; return; }

    local wall rss p1 p2 n_parts unique
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")
    n_parts=$(se_val "n_parts" "$se")
    unique=$(se_val "unique_kmers" "$se")

    printf "  %-10s %-8s n=%-3d t=%-3d ram=%-4sGB  wall=%8.3fs  RSS=%s MB\n" \
        "$tool" "$ds" "$n_files" "$threads" "$ram_gb" "$wall" "$rss"
    echo "$ds,$n_files,$threads,$ram_gb,$tool,$wall,$rss,$p1,$p2,$n_parts,$unique" >> "$CSV"
}

# =============================================================================
# Main loop: each experiment (ram_gb + tool subset) × both datasets
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"
for exp in "${EXPERIMENTS[@]}"; do
    ram_gb="${exp%%:*}"
    tools="${exp#*:}"
    echo "── ram=${ram_gb}GB  tools=${tools}  t=${THREADS} ──────────────────────"
    for tool in $tools; do
        run_one "$tool" gallus "$GALLUS_FOF" "$GALLUS_N" "$THREADS" "$ram_gb"
        run_one "$tool" human3 "$HUMAN3_FOF" "$HUMAN3_N" "$THREADS" "$ram_gb"
    done
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"