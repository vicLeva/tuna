#!/usr/bin/env bash
# scripts/bench_scaling.sh — n-files scaling sweep, count-only, 3 tools
#
# Sweeps the number of input files for tuna, KMC, and FastK on E. coli and
# Human assembly datasets.  Tool order: one tool processes all (dataset × N)
# combinations before the next tool starts — so partial runs give complete
# curves per tool.
#
# Count-only for all 3 tools — no disk output, wall-time reflects counting
# alone (each tool's own native mechanism, verified against source/docs):
#   tuna  : -co (skip output writing after counting)
#   KMC   : -w (without output) — binary db never written; kmc_dump step
#           removed entirely (nothing to dump)
#   FastK : no -t/-p — DO_TABLE/DO_PROFILE stay 0, so Merge_Tables/
#           Merge_Profiles never run (core Sorting() still happens
#           unconditionally); Tabex dump step removed entirely (no .ktab
#           ever gets written, so there is nothing to dump)
#
# Phase breakdown:
#   tuna : own "phase1:"/"phase2:" stderr lines
#   KMC  : "1st stage:"/"2nd stage:" stderr lines — printed unconditionally
#          by kmc_CLI/kmc.cpp's print_summary() (undocumented in -h usage
#          text, but always emitted; verified directly against KMC's own
#          CLI source and empirically against the local binary)
#   FastK: no equivalent phase breakdown found; left as "na"
#
# Usage:
#   bash scripts/bench_scaling.sh
#   THREADS=16 bash scripts/bench_scaling.sh
#
# Output: $RESULTS/bench_scaling.csv
#   dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s
#
# Notes:
#   FastK does not support @fof — files are passed as individual arguments;
#   .fna files need .fa symlinks (created once at setup in $WORKDIR/fastk_links/)

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================

: "${TUNA:=""}"      # /path/to/tuna
: "${KMC:=""}"       # /path/to/kmc
: "${FASTK:=""}"     # /path/to/FastK

WORKDIR="/WORKS/vlevallois/expes_tuna"
ECOLI_FOF="/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list"
HUMAN_FOF="/WORKS/vlevallois/data/dataset_genome_human/fof.list"

K=31
M=21
RAM_GB=256
THREADS=${THREADS:-8}

ECOLI_NS=(1 2 3 5 10 20 50 100 200 500 1000 1500 2000 2500 3000 3500)
HUMAN_NS=(1 2 3 4 5 6 7 8 9 10 15 20 25 30 60)

TOOLS=(tuna kmc fastk)

# =============================================================================
# Sanity checks
# =============================================================================

err=0
for var in TUNA KMC FASTK; do
    [[ -z "${!var}" ]] && { echo "[error] $var is not set"; err=1; }
done
for var in TUNA KMC FASTK; do
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] not executable: ${!var}"; err=1; }
done
for var in ECOLI_FOF HUMAN_FOF; do
    [[ -n "${!var}" && ! -f "${!var}" ]] && { echo "[error] not found: ${!var}"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

ECOLI_TOTAL=$(wc -l < "$ECOLI_FOF")
HUMAN_TOTAL=$(wc -l < "$HUMAN_FOF")
echo "[info] E. coli fof: $ECOLI_TOTAL files"
echo "[info] Human   fof: $HUMAN_TOTAL files"

ECOLI_MAX=${ECOLI_NS[-1]}
HUMAN_MAX=${HUMAN_NS[-1]}
[[ "$ECOLI_MAX" -gt "$ECOLI_TOTAL" ]] && \
    echo "[warn] largest E. coli N ($ECOLI_MAX) > fof size ($ECOLI_TOTAL); those points will be skipped"
[[ "$HUMAN_MAX" -gt "$HUMAN_TOTAL" ]] && \
    echo "[warn] largest Human N ($HUMAN_MAX) > fof size ($HUMAN_TOTAL); those points will be skipped"

# =============================================================================
# Setup
# =============================================================================

RESULTS="$WORKDIR/bench_scaling_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS" "$WORKDIR/fastk_links"

CSV="$RESULTS/bench_scaling.csv"
echo "dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s" > "$CSV"
echo "[bench] Results: $CSV"
echo "[bench] k=$K  m=$M  threads=$THREADS  ram=${RAM_GB}GB  (count-only)"

# FastK requires .fa / .fa.gz extensions; .fna files need symlinks.
# Build FastK-compatible fofs once upfront for all files we will use.

fastk_link() {
    local f="$1"
    local base; base=$(basename "$f")
    case "$base" in
        *.fna.gz) echo "$WORKDIR/fastk_links/${base%.fna.gz}.fa.gz" ;;
        *.fna)    echo "$WORKDIR/fastk_links/${base%.fna}.fa"       ;;
        *)        echo "$WORKDIR/fastk_links/$base"                  ;;
    esac
}

echo "[setup] Creating FastK symlinks..."
ECOLI_FASTK_FOF="$RESULTS/fastk_ecoli_full.fof"
HUMAN_FASTK_FOF="$RESULTS/fastk_human_full.fof"
> "$ECOLI_FASTK_FOF"
> "$HUMAN_FASTK_FOF"

while IFS= read -r f; do
    lnk=$(fastk_link "$f")
    ln -sf "$f" "$lnk"
    echo "$lnk" >> "$ECOLI_FASTK_FOF"
done < <(head -n "$ECOLI_MAX" "$ECOLI_FOF" 2>/dev/null || head -n "$ECOLI_TOTAL" "$ECOLI_FOF")

while IFS= read -r f; do
    lnk=$(fastk_link "$f")
    ln -sf "$f" "$lnk"
    echo "$lnk" >> "$HUMAN_FASTK_FOF"
done < <(head -n "$HUMAN_MAX" "$HUMAN_FOF" 2>/dev/null || head -n "$HUMAN_TOTAL" "$HUMAN_FOF")

echo "[setup] Done."

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
# Run functions — count-only, single timed call per tool, no dump step
# =============================================================================

run_tuna() {
    local ds="$1" n="$2" fof="$3"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/subfof_${ds}_${n}.list"
    head -n "$n" "$fof" > "$subfof"

    local work="$WORKDIR/tuna_scaling_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_tuna.timefile"
    local se="$RESULTS/${ds}_n${n}_tuna.stderr"
    mkdir -p "$work"

    # -co: skip output writing after counting
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -hp -co \
        -ram "$RAM_GB" \
        -w "$work/" "@$subfof" /dev/null \
        > /dev/null 2>"$se" \
    || { echo "  [FAIL] tuna $ds n=$n"; rm -rf "$work" "$subfof"; return; }
    rm -rf "$work" "$subfof"

    local wall rss p1 p2
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")

    printf "    tuna   n=%-5d  wall=%7ss  p1=%ss  p2=%ss  RSS=%sMB\n" \
        "$n" "$wall" "$p1" "$p2" "$rss"
    echo "$ds,$n,tuna,$wall,$rss,$p1,$p2" >> "$CSV"
}

run_kmc() {
    local ds="$1" n="$2" fof="$3"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/subfof_${ds}_${n}.list"
    head -n "$n" "$fof" > "$subfof"

    local db="$WORKDIR/kmc_scaling_${ds}_n${n}"
    local tmp="$WORKDIR/kmc_scaling_tmp_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_kmc.timefile"
    local se="$RESULTS/${ds}_n${n}_kmc.stderr"
    mkdir -p "$tmp"

    # -w: without output — binary db is never written (true count-only)
    # Note: this cluster's KMC build prints its 1st/2nd-stage summary to
    # stdout, not stderr — capture both streams into $se so se_val can find it.
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k"$K" -m"$RAM_GB" -ci1 -cs4294967295 -fm -hp -t"$THREADS" -w \
        "@$subfof" "$db" "$tmp" \
        > "$se" 2>&1 \
    || { echo "  [FAIL] kmc $ds n=$n"; rm -rf "$tmp" "$subfof"; return; }
    rm -rf "$tmp" "$subfof"

    local wall rss p1 p2
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "1st stage" "$se")
    p2=$(se_val "2nd stage" "$se")

    printf "    kmc    n=%-5d  wall=%7ss  p1=%ss  p2=%ss  RSS=%sMB\n" \
        "$n" "$wall" "$p1" "$p2" "$rss"
    echo "$ds,$n,kmc,$wall,$rss,$p1,$p2" >> "$CSV"
}

run_fastk() {
    local ds="$1" n="$2" fof="$3" fastk_fof="$4"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/fastk_subfof_${ds}_${n}.list"
    head -n "$n" "$fastk_fof" > "$subfof"

    local db="$WORKDIR/fastk_scaling_${ds}_n${n}"
    local tmp="$WORKDIR/fastk_scaling_tmp_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_fastk.timefile"
    local se="$RESULTS/${ds}_n${n}_fastk.stderr"
    mkdir -p "$tmp"

    # No -t/-p: DO_TABLE/DO_PROFILE stay 0, so no .ktab ever gets written —
    # true count-only (core Sorting() still happens unconditionally either
    # way); no Tabex dump step since there is nothing to dump.
    mapfile -t fastk_files < "$subfof"
    /usr/bin/time -v -o "$tf" \
        "$FASTK" -k"$K" -T"$THREADS" -M"$RAM_GB" \
        -N"$db" -P"$tmp" "${fastk_files[@]}" \
        > /dev/null 2>"$se" \
    || { echo "  [FAIL] FastK $ds n=$n"; rm -rf "$tmp" "$subfof"; return; }
    rm -rf "$tmp"
    rm -f "${db}"* "$subfof" 2>/dev/null || true

    local wall rss
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")

    printf "    fastk  n=%-5d  wall=%7ss  RSS=%sMB\n" "$n" "$wall" "$rss"
    echo "$ds,$n,fastk,$wall,$rss,na,na" >> "$CSV"
}

# =============================================================================
# Main loop: one full tool before moving to the next
# =============================================================================

echo ""
echo "[bench] Starting — $(date)"
echo ""

for tool in "${TOOLS[@]}"; do
    echo "════ Tool: $tool ════════════════════════════════════════════"

    echo "  ── ecoli ──"
    for n in "${ECOLI_NS[@]}"; do
        case "$tool" in
            tuna)  run_tuna  ecoli "$n" "$ECOLI_FOF" ;;
            kmc)   run_kmc   ecoli "$n" "$ECOLI_FOF" ;;
            fastk) run_fastk ecoli "$n" "$ECOLI_FOF" "$ECOLI_FASTK_FOF" ;;
        esac
    done

    echo "  ── human ──"
    for n in "${HUMAN_NS[@]}"; do
        case "$tool" in
            tuna)  run_tuna  human "$n" "$HUMAN_FOF" ;;
            kmc)   run_kmc   human "$n" "$HUMAN_FOF" ;;
            fastk) run_fastk human "$n" "$HUMAN_FOF" "$HUMAN_FASTK_FOF" ;;
        esac
    done

    echo ""
done

echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
echo ""
column -t -s, "$CSV" | head -20
echo "..."
