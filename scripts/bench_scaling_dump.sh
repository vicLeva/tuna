#!/usr/bin/env bash
# scripts/bench_scaling_dump.sh — n-files scaling sweep, full dump, 3 tools
#
# Bis of bench_scaling.sh with real output written to disk instead of
# count-only, to measure the actual dump cost (vs. new_runs_v2's count-only
# baseline).
#
# Sweeps the number of input files for tuna, KMC, and FastK on E. coli and
# Human assembly datasets. Tool order: one tool processes all (dataset × N)
# combinations before the next tool starts.
#
# Full dump for all 3 tools — output files are deleted right after each
# run's timing is captured (not kept — full dumps can be very large):
#   tuna  : writes its TSV (no -co). Dump cost is interleaved per-partition
#           into phase2 and can't be isolated directly here — compare
#           phase2_s against new_runs_v2's count-only phase2_s to estimate it.
#   KMC   : writes its binary db (no -w), then a separately-timed kmc_dump
#           step (dump_s) to ASCII. KMC's own binary-write cost is folded
#           into its stage times same as tuna — compare 1st/2nd stage here
#           against new_runs_v2 to estimate it. dump_s is the extra ASCII
#           dump step only.
#   FastK : -t1 to produce the .ktab table, then Tabex -A to dump ASCII,
#           timed together as one combined wall time (matches
#           bench_all_tools.sh's convention). Expect frequent failures —
#           Tabex crashes are common; failures are logged and skipped, not
#           fatal to the sweep. Also still subject to the known
#           Scan_All_Input N=threads+1 hang (root-caused 2026-07-21,
#           confirmed independent of -t/-p) — kept under the same timeout.
#
# Phase breakdown:
#   tuna : own "phase1:"/"phase2:" stderr lines
#   KMC  : "1st stage:"/"2nd stage:" stderr lines
#   FastK: no equivalent phase breakdown; left as "na"
#
# Usage:
#   bash scripts/bench_scaling_dump.sh
#   THREADS=16 bash scripts/bench_scaling_dump.sh
#
# Output: $RESULTS/bench_scaling_dump.csv
#   dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,dump_s
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
: "${KMC_DUMP:=""}"  # /path/to/kmc_dump
: "${FASTK:=""}"     # /path/to/FastK
: "${TABEX:=""}"     # /path/to/Tabex

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
for var in TUNA KMC KMC_DUMP FASTK TABEX; do
    [[ -z "${!var}" ]] && { echo "[error] $var is not set"; err=1; }
done
for var in TUNA KMC KMC_DUMP FASTK TABEX; do
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

RESULTS="$WORKDIR/bench_scaling_dump_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS" "$WORKDIR/fastk_links"

CSV="$RESULTS/bench_scaling_dump.csv"
echo "dataset,n_files,tool,wall_s,rss_mb,phase1_s,phase2_s,dump_s" > "$CSV"
echo "[bench] Results: $CSV"
echo "[bench] k=$K  m=$M  threads=$THREADS  ram=${RAM_GB}GB  (full dump)"

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
# Run functions — full dump
# =============================================================================

run_tuna() {
    local ds="$1" n="$2" fof="$3"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/subfof_${ds}_${n}.list"
    head -n "$n" "$fof" > "$subfof"

    local work="$WORKDIR/tuna_scaling_dump_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_tuna.timefile"
    local se="$RESULTS/${ds}_n${n}_tuna.stderr"
    mkdir -p "$work"

    # Full dump: real TSV output instead of /dev/null, no -co.
    /usr/bin/time -v -o "$tf" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" -hp \
        -ram "$RAM_GB" \
        -w "$work/" "@$subfof" "$work/out.tsv" \
        > /dev/null 2>"$se"
    local rc=$?
    rm -f "$work/out.tsv"
    if [[ "$rc" -ne 0 ]]; then echo "  [FAIL] tuna $ds n=$n"; rm -rf "$work" "$subfof"; return; fi
    rm -rf "$work" "$subfof"

    local wall rss p1 p2
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "phase1" "$se")
    p2=$(se_val "phase2" "$se")

    printf "    tuna   n=%-5d  wall=%7ss  p1=%ss  p2=%ss  RSS=%sMB\n" \
        "$n" "$wall" "$p1" "$p2" "$rss"
    echo "$ds,$n,tuna,$wall,$rss,$p1,$p2,na" >> "$CSV"
}

run_kmc() {
    local ds="$1" n="$2" fof="$3"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/subfof_${ds}_${n}.list"
    head -n "$n" "$fof" > "$subfof"

    local db="$WORKDIR/kmc_scaling_dump_${ds}_n${n}"
    local tmp="$WORKDIR/kmc_scaling_dump_tmp_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_kmc.timefile"
    local se="$RESULTS/${ds}_n${n}_kmc.stderr"
    local dt="$RESULTS/${ds}_n${n}_kmc.dumptime"
    mkdir -p "$tmp"

    # Full dump: real binary db (no -w). This cluster's KMC build prints its
    # 1st/2nd-stage summary to stdout, not stderr — capture both streams.
    /usr/bin/time -v -o "$tf" \
        "$KMC" -k"$K" -m"$RAM_GB" -ci1 -cs4294967295 -fm -hp -t"$THREADS" \
        "@$subfof" "$db" "$tmp" \
        > "$se" 2>&1
    local rc=$?
    rm -rf "$tmp"
    if [[ "$rc" -ne 0 ]]; then echo "  [FAIL] kmc $ds n=$n"; rm -rf "$subfof" "${db}"* 2>/dev/null; return; fi

    /usr/bin/time -v -o "$dt" \
        "$KMC_DUMP" "$db" "$db.tsv" \
        > /dev/null 2>>"$se"
    rc=$?
    rm -f "$subfof" "$db.tsv" "${db}"* 2>/dev/null
    if [[ "$rc" -ne 0 ]]; then echo "  [FAIL] kmc_dump $ds n=$n"; return; fi

    local wall rss p1 p2 dump
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")
    p1=$(se_val "1st stage" "$se")
    p2=$(se_val "2nd stage" "$se")
    dump=$(wall_from_file "$dt")

    printf "    kmc    n=%-5d  wall=%7ss  p1=%ss  p2=%ss  dump=%ss  RSS=%sMB\n" \
        "$n" "$wall" "$p1" "$p2" "$dump" "$rss"
    echo "$ds,$n,kmc,$wall,$rss,$p1,$p2,$dump" >> "$CSV"
}

run_fastk() {
    local ds="$1" n="$2" fof="$3" fastk_fof="$4"
    local total; total=$(wc -l < "$fof")
    if [[ "$n" -gt "$total" ]]; then echo "  [skip] $ds n=$n > fof ($total)"; return; fi

    local subfof="$RESULTS/fastk_subfof_${ds}_${n}.list"
    head -n "$n" "$fastk_fof" > "$subfof"

    local db="$WORKDIR/fastk_scaling_dump_${ds}_n${n}"
    local tmp="$WORKDIR/fastk_scaling_dump_tmp_${ds}_n${n}"
    local tf="$RESULTS/${ds}_n${n}_fastk.timefile"
    local se="$RESULTS/${ds}_n${n}_fastk.stderr"
    mkdir -p "$tmp"

    # -t1: produce the .ktab table (count >= 1), then Tabex -A dumps ASCII.
    # Both timed together as one combined wall time (matches
    # bench_all_tools.sh's convention). Expect frequent failures here —
    # Tabex crashes on a meaningful fraction of runs.
    #
    # Same known FastK bug as bench_scaling.sh (Scan_All_Input hangs when
    # input file count == thread count + 1, root-caused 2026-07-21,
    # confirmed independent of -t/-p) — kept under the same timeout so one
    # bad (n, THREADS) pair can't stall the whole sweep.
    mapfile -t fastk_files < "$subfof"
    timeout 900 /usr/bin/time -v -o "$tf" bash -c "
        set -e
        \"$FASTK\" -k${K} -t1 -T${THREADS} -M${RAM_GB} \
            -N\"$db\" -P\"$tmp\" ${fastk_files[*]@Q} \
        && \"$TABEX\" -A \"$db\" > \"$db.tsv\"
    " > /dev/null 2>"$se"
    local rc=$?
    rm -rf "$tmp" "$db.tsv" "${db}"* 2>/dev/null
    if [[ "$rc" -eq 124 ]]; then
        echo "  [TIMEOUT] FastK $ds n=$n (known N=threads+1 hang, see memory)"
        echo "$ds,$n,fastk,timeout,,na,na,na" >> "$CSV"
        rm -f "$subfof"
        return
    elif [[ "$rc" -ne 0 ]]; then
        echo "  [FAIL] FastK $ds n=$n (Tabex crash or FastK error)"
        rm -f "$subfof"
        return
    fi
    rm -f "$subfof"

    local wall rss
    wall=$(wall_from_file "$tf")
    rss=$(rss_mb "$tf")

    printf "    fastk  n=%-5d  wall=%7ss  RSS=%sMB\n" "$n" "$wall" "$rss"
    echo "$ds,$n,fastk,$wall,$rss,na,na,na" >> "$CSV"
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
