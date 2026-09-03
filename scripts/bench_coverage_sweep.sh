#!/usr/bin/env bash
# bench_coverage_sweep.sh - coverage sweep, tuna vs KMC
#
# This script produced benchmark/results/coverage_sweep_kh/ on 2026-08-28 at
# m=17, and coverage_sweep_ud/ afterwards at m=21. The default is now m=21 so
# the coverage curve is comparable with the rest of the suite; the older CSV
# records m=17 in its own m column, so the two remain distinguishable.
#
# It still does not share bench_common.sh, and its default ROOT still points at
# expes_tuna/coverage_sweep rather than expes_paper: the rest of its design is
# specific to this experiment and is left as it was.
#
# Source: SRR622461_1.fastq.gz, human reads, measured coverage ~5x.
#
# Building the levels
#   Below 5x these are genuine subsamples: level Nx takes the first N/5 of the
#   reads, so k-mer diversity really does grow the way it would with shallower
#   sequencing, and 5x reproduces the original file exactly.
#   Above 5x there are no further reads to draw on, so levels are duplicated:
#   Nx = (N div 5) copies of the full file + the first (N mod 5)/5 of it.
#   gzip streams concatenate natively, so nothing is recompressed beyond the
#   four fractional blocks, which are built once and reused.
#
#   Keep in mind: above 5x the number of DISTINCT k-mers stops growing, only
#   the counts rise. These levels measure how the tools scale with input
#   volume, not with true sequencing depth.
#
# Order of work: one full pass per tool, rather than interleaving.
#   pass 1  tuna bin    (KFF output)     over every level
#   pass 2  tuna ascii  (TSV output)     over every level
#   pass 3  kmc  bin    (binary DB, then kmc_dump timed separately)
#   Running a whole pass per tool keeps conditions uniform within a tool: no
#   tool ever lands right after another has freed tens of GB of page cache.
#
# Level files inside a pass are built incrementally. Going from Nx to (N+5)x
# is a rename plus one appended copy of the source, so each level is written
# once instead of being rebuilt from scratch: ~55 GB of writes per pass rather
# than several hundred. The file for level N is consumed by level N+5, so at
# most one level file exists at a time. Each pass rebuilds the chain from the
# first level it actually needs.
#
# The ASCII price is read afterwards: KMC = dump_s, tuna = ascii - bin.
#
# TSV outputs are by far the biggest (~40 bytes per distinct k-mer, against ~5
# for binary). Before each ASCII step - kmc_dump as well as tuna ascii - the
# remaining space is checked: if it is not enough the step is skipped (status
# "skipped_space") and the level carries on instead of failing. Binary runs
# always go ahead.
#
# Stops as soon as a level would exceed STOP_BYTES (500 GB by default).
# Resumes at measurement granularity: anything already in the CSV is skipped.
#
# Statistics (when KMC_TOOLS is set)
#   The count histogram is taken from the KMC binary DB, which exists at that
#   point anyway: ~5 bytes per distinct k-mer against ~40 for the TSV. Past
#   twenty-odd x the TSV runs to terabytes and a pass over it would cost hours;
#   the histogram stays in the megabytes. Everything else, skew included,
#   follows from it. Statistics are therefore produced during the KMC pass.
#
# Output: $ROOT/coverage_sweep.csv
#   coverage,m,file_bytes,tool,mode,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers,total_kmers,status

set -uo pipefail
export LC_ALL=C LANG=C

: "${TUNA:=}"
: "${KMC:=}"
: "${KMC_DUMP:=}"
: "${KMC_TOOLS:=}"        # optional: enables per-level statistics
: "${ROOT:=/WORKS/vlevallois/expes_tuna/coverage_sweep}"

: "${SRC:=/WORKS/vlevallois/data/data_sequencing/SRR622461_1.fastq.gz}"
SRC_COV=5                                  # coverage of the source file
K=31
M=${M:-21}          # m=21 everywhere in the suite, so the coverage curve is
                    # comparable with the other experiments. m=17 was measured
                    # -3.8% time and -15% RSS on these reads and was the
                    # original default here, which is why the first run of this
                    # sweep (2026-08-28) carries m=17 in its CSV.
THREADS=${THREADS:-8}
RAM_GB=${RAM_GB:-256}
STOP_BYTES=${STOP_BYTES:-$((500 * 1024 * 1024 * 1024))}
# Default levels: 1..10, then in steps of 5 up to 100.
# Override without editing the script:  COV_LIST="1 5 10 50 100" bash ...
if [[ -n "${COV_LIST:-}" ]]; then
    read -ra COVERAGES <<< "$COV_LIST"
else
    COVERAGES=(1 2 3 4 5 6 7 8 9 10 $(seq 15 5 100))
fi
# Tool passes, in order. Each entry is "<tool> <mode>".
PASSES=("tuna bin" "tuna ascii" "kmc bin")

WORK="$ROOT/work"; AUX="$ROOT/aux"; DATA="$ROOT/data"; STATS="$ROOT/stats"
CSV="$ROOT/coverage_sweep.csv"
SCSV="$ROOT/coverage_stats.csv"
CHUNK_DIR="$DATA"          # fractional blocks chunk_{1..4}_5.fastq.gz

err=0
for v in TUNA KMC KMC_DUMP; do
    [[ -z "${!v}" ]]              && { echo "[error] $v is not set"; err=1; }
    [[ -n "${!v}" && ! -x "${!v}" ]] && { echo "[error] not executable: ${!v}"; err=1; }
done
[[ -f "$SRC" ]] || { echo "[error] not found: $SRC"; err=1; }
[[ $err -eq 1 ]] && exit 1

mkdir -p "$WORK" "$AUX" "$DATA" "$STATS"
[[ -f "$CSV" ]] || echo "coverage,m,file_bytes,tool,mode,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers,total_kmers,status" > "$CSV"
[[ -f "$SCSV" ]] || echo "coverage,m,file_bytes,distinct_kmers,total_kmers,singletons,pct_singletons,mean_count,median_count,p99_count,max_count,cov_peak,top1pct_mass_pct,gini" > "$SCSV"

echo "[bench] coverage sweep - k=$K m=$M t=$THREADS ram=${RAM_GB}GB"
echo "[bench] source : $SRC ($(du -h "$SRC" | cut -f1), ~${SRC_COV}x)"
echo "[bench] tuna   : $TUNA  ($(stat -c %y "$TUNA" 2>/dev/null | cut -d. -f1))"
echo "[bench] kmc    : $KMC"
echo "[bench] output : $CSV"

# -- helpers ----------------------------------------------------------------
wall_of() { awk -F: '/Elapsed \(wall clock\)/{t=$0} END{n=split(t,a," ");split(a[n],b,":");
            if(length(b)==3) printf "%.3f",b[1]*3600+b[2]*60+b[3]; else printf "%.3f",b[1]*60+b[2]}' "$1"; }
rss_of()  { awk '/Maximum resident/{printf "%.0f",$NF/1024}' "$1"; }
se_val()  { grep -m1 "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2;gsub(/^ +/,"",v);sub(/s$/,"",v);print v}'; }
kmc_num() { grep -m1 "$1" "$2" 2>/dev/null | awk -F: '{gsub(/^ +/,"",$2);print $2}'; }
free_bytes() { df -B1 --output=avail "$1" | tail -1; }
# TSV ~ (k + separator + count + newline) bytes per distinct k-mer
tsv_bytes_for() { echo $(( $1 * (K + 12) )); }
room_for() { local a; a=$(free_bytes "$WORK"); (( a > $1 * 115 / 100 )); }
row()     { echo "$1" >> "$CSV"; }
# Resume at measurement granularity: an interrupted pass picks up exactly
# where it stopped instead of skipping levels it only partly covered.
have_run() { awk -F, -v c="$1" -v t="$2" -v m="$3" \
    'NR>1 && $1==c && $4==t && $5==m {f=1} END{exit !f}' "$CSV"; }
have_stats() { awk -F, -v c="$1" 'NR>1 && $1==c {f=1} END{exit !f}' "$SCSV"; }
# Distinct-k-mer count already recorded for this level, whatever produced it.
# Used to size the TSV space check when the ASCII pass runs on its own.
distinct_for() { awk -F, -v c="$1" 'NR>1 && $1==c && $11!="" {v=$11} END{print v}' "$CSV"; }

# -- per-level statistics, from the count histogram -------------------------
# hist: one "<count> <number of distinct k-mers with that count>" per line.
kmer_stats() {
    local cov="$1" fb="$2" db="$3"
    [[ -z "$KMC_TOOLS" ]] && return
    have_stats "$cov" && return
    [[ -s "${db}.kmc_suf" ]] || { echo "    [stats] no KMC DB, skipped"; return; }
    local hist="$STATS/hist_cov${cov}x.tsv"
    "$KMC_TOOLS" transform "$db" histogram "$hist" -ci1 -cx1000000 >/dev/null 2>&1 \
        || { echo "    [stats] histogram unavailable"; return; }
    awk -v cov="$cov" -v mm="$M" -v fb="$fb" '
        { c=$1+0; n=$2+0; if(n==0) next
          cnt[++m]=c; num[m]=n; d+=n; t+=c*n; if(c==1) s1=n; if(c>mx) mx=c
          if(c>1 && n>peakn){peakn=n; peak=c} }
        END{
          if(d==0){print cov","mm","fb",,,,,,,,,,"; exit}
          # medians / quantiles over distinct k-mers, counts ascending
          need50=d*0.5; need99=d*0.99; need1=d*0.01
          cum=0; for(i=1;i<=m;i++){ cum+=num[i]
            if(!med && cum>=need50) med=cnt[i]
            if(!p99 && cum>=need99) p99=cnt[i] }
          # mass held by the top 1% most abundant k-mers -> skew
          cum=0; mass=0; for(i=m;i>=1;i--){ take=num[i]
            if(cum+take>need1) take=need1-cum
            mass+=take*cnt[i]; cum+=take; if(cum>=need1) break }
          # Gini over the count distribution
          cum=0; area=0; for(i=1;i<=m;i++){ area+=num[i]*(2*cum+num[i])*cnt[i]; cum+=num[i] }
          gini=(d>1&&t>0)? area/(d*t)-1 : 0
          printf "%s,%s,%s,%d,%d,%d,%.3f,%.3f,%d,%d,%d,%d,%.2f,%.4f\n",
                 cov,mm,fb,d,t,s1,100*s1/d,t/d,med,p99,mx,peak,100*mass/t,gini
        }' "$hist" >> "$SCSV"
    echo "    [stats] $(tail -1 "$SCSV" | awk -F, '{printf "distinct=%s  singletons=%.1f%%  peak=%sx  max=%s  top1%%=%s%%  gini=%s",$4,$7,$12,$11,$13,$14}')"
}

# -- fractional blocks 1/5 .. 4/5 (built once) ------------------------------
LINES_FILE="$DATA/.total_lines"
if [[ ! -s "$LINES_FILE" ]]; then
    echo "[setup] counting reads in the source (one decompression pass)..."
    zcat "$SRC" | wc -l > "$LINES_FILE"
fi
TOTAL_LINES=$(cat "$LINES_FILE")
LINES_1X=$(( (TOTAL_LINES / SRC_COV) / 4 * 4 ))     # aligned on 4 lines = 1 read
echo "[setup] $TOTAL_LINES lines -> $LINES_1X lines per fifth ($((LINES_1X/4)) reads)"

ZIP=$(command -v pigz >/dev/null && echo "pigz -1 -p $THREADS" || echo "gzip -1")
for r in 1 2 3 4; do
    ck="$CHUNK_DIR/chunk_${r}_5.fastq.gz"
    [[ -s "$ck" ]] && continue
    echo "[setup] building block ${r}/5 with: $ZIP"
    # head closes the pipe: zcat takes a SIGPIPE and the pipeline returns
    # non-zero. That is expected, so pipefail is turned off here and success
    # is judged on the resulting file.
    set +o pipefail
    zcat "$SRC" | head -n $(( LINES_1X * r )) | $ZIP > "$ck.tmp"
    set -o pipefail
    [[ -s "$ck.tmp" ]] || { echo "[error] block ${r}/5 is empty"; rm -f "$ck.tmp"; exit 1; }
    [[ $(zcat "$ck.tmp" 2>/dev/null | head -n 1 | cut -c1) == "@" ]] \
        || { echo "[error] block ${r}/5 is not FASTQ"; rm -f "$ck.tmp"; exit 1; }
    mv "$ck.tmp" "$ck"
done
SRC_BYTES=$(stat -c%s "$SRC")

CHUNK_BYTES=$(( SRC_BYTES / SRC_COV ))
echo "[setup] 1x ~ $(numfmt --to=iec "$CHUNK_BYTES")  (=> 100x ~ $(numfmt --to=iec $((CHUNK_BYTES*100))))"

# -- level construction ------------------------------------------------------
# A level is q copies of the source plus the first r fifths, with q = cov/5 and
# r = cov%5. Going from a level with r == 0 to a larger one with r == 0 is
# therefore just "append q'-q copies", which is done by renaming the previous
# file and appending to it. Every other transition is rebuilt from scratch;
# those are the small levels (below 10x, at most 4.5 GB), so it costs nothing.
#
# args: target_cov previous_cov ("" if none)   -> leaves $DATA/cov_<target>x.fastq.gz
build_level() {
    local cov="$1" prev="${2:-}"
    local f="$DATA/cov_${cov}x.fastq.gz"
    local pf="" i
    [[ -n "$prev" ]] && pf="$DATA/cov_${prev}x.fastq.gz"
    # Already there (leftover from an interrupted run): reuse it, but still
    # drop the predecessor, otherwise level files pile up and the "one file at
    # a time" space budget no longer holds.
    if [[ -s "$f" ]]; then
        [[ -n "$pf" && -s "$pf" && "$pf" != "$f" ]] && rm -f "$pf"
        return 0
    fi

    local q=$(( cov / SRC_COV )) r=$(( cov % SRC_COV ))
    if [[ -n "$prev" && -s "$pf" ]] && (( r == 0 && prev % SRC_COV == 0 && prev < cov )); then
        local pq=$(( prev / SRC_COV ))
        echo "    building ${cov}x from ${prev}x (+$(( q - pq )) copies)"
        mv "$pf" "$f" || return 1
        for ((i=pq;i<q;i++)); do cat "$SRC" >> "$f" || return 1; done
    else
        [[ -n "$pf" && -s "$pf" ]] && rm -f "$pf"     # chain broken, drop it
        echo "    building ${cov}x from scratch"
        : > "$f"
        for ((i=0;i<q;i++)); do cat "$SRC" >> "$f" || return 1; done
        (( r > 0 )) && { cat "$CHUNK_DIR/chunk_${r}_5.fastq.gz" >> "$f" || return 1; }
    fi
    [[ -s "$f" ]] || return 1
    return 0
}

# -- one measurement ---------------------------------------------------------
# args: cov tool mode
run_one() {
    local cov="$1" tool="$2" mode="$3"
    local f="$DATA/cov_${cov}x.fastq.gz"
    local fb; fb=$(stat -c%s "$f" 2>/dev/null || echo 0)
    local tag="cov${cov}_${tool}_${mode}"
    local tf="$AUX/$tag.time" se="$AUX/$tag.stderr" dt="$AUX/$tag.dumptime"
    local wall rss p1 p2 dump uniq tot rc
    # KMC leaves a huge RSS that evicts the page cache: without this read the
    # first tuna run would read cold and the second warm, skewing the bin/ascii
    # gap. So warm up before every measurement.
    cat "$f" > /dev/null 2>&1

    case "$tool:$mode" in
      kmc:bin)
        mkdir -p "$WORK/tmp"
        /usr/bin/time -v -o "$tf" "$KMC" -k$K -ci1 -cs4294967295 -fq -m${RAM_GB} -hp -t${THREADS} \
            "$f" "$WORK/out" "$WORK/tmp" > "$se" 2>&1
        rc=$?; rm -rf "$WORK/tmp"
        if [[ $rc -ne 0 ]]; then row "$cov,na,$fb,kmc,bin,,,,,,,,fail"; echo "    [FAIL] kmc bin"; return; fi
        uniq=$(kmc_num "No. of unique k-mers" "$se")
        local want; want=$(tsv_bytes_for "${uniq:-0}")
        if [[ -n "$uniq" ]] && ! room_for "$want"; then
            dump="skipped"
            echo "    [skip] kmc_dump: would need ~$(numfmt --to=iec "$want"), free $(numfmt --to=iec "$(free_bytes "$WORK")")"
        else
            /usr/bin/time -v -o "$dt" "$KMC_DUMP" "$WORK/out" "$WORK/out.tsv" >/dev/null 2>>"$se"
            rc=$?; dump=$([[ $rc -eq 0 ]] && wall_of "$dt" || echo "")
            rm -f "$WORK/out.tsv"
        fi
        kmer_stats "$cov" "$fb" "$WORK/out"     # before dropping the DB
        rm -f "$WORK/out.kmc_pre" "$WORK/out.kmc_suf"
        wall=$(wall_of "$tf"); rss=$(rss_of "$tf")
        p1=$(se_val "1st stage" "$se"); p2=$(se_val "2nd stage" "$se")
        tot=$(kmc_num "Total no. of k-mers" "$se")
        row "$cov,na,$fb,kmc,bin,$wall,$rss,$p1,$p2,$dump,$uniq,$tot,ok"
        printf "    [%-4s %-5s] wall=%9ss  p1=%9ss  p2=%9ss  dump=%9ss  RSS=%8sMB\n" kmc bin "$wall" "$p1" "$p2" "$dump" "$rss"
        ;;
      tuna:*)
        if [[ $mode == ascii ]]; then
            local est; est=$(distinct_for "$cov")
            if [[ -n "$est" ]]; then
                local want; want=$(tsv_bytes_for "$est")
                if ! room_for "$want"; then
                    row "$cov,$M,$fb,tuna,ascii,,,,,,,,skipped_space"
                    echo "    [skip] tuna ascii: would need ~$(numfmt --to=iec "$want"), free $(numfmt --to=iec "$(free_bytes "$WORK")")"
                    return
                fi
            fi
        fi
        local out; [[ $mode == bin ]] && out="$WORK/out.kff" || out="$WORK/out.tsv"
        /usr/bin/time -v -o "$tf" "$TUNA" -k $K -m $M -t $THREADS -ram $RAM_GB -hp \
            -w "$WORK/" "$f" "$out" >/dev/null 2>"$se"
        rc=$?; rm -f "$out"
        if [[ $rc -ne 0 ]]; then row "$cov,$M,$fb,tuna,$mode,,,,,,,,fail"; echo "    [FAIL] tuna $mode"; return; fi
        wall=$(wall_of "$tf"); rss=$(rss_of "$tf")
        p1=$(se_val phase1 "$se"); p2=$(se_val phase2 "$se")
        uniq=$(se_val unique_kmers "$se"); tot=$(se_val total_kmers "$se")
        row "$cov,$M,$fb,tuna,$mode,$wall,$rss,$p1,$p2,,$uniq,$tot,ok"
        printf "    [%-4s %-5s] wall=%9ss  p1=%9ss  p2=%9ss  RSS=%8sMB  uniq=%s\n" tuna "$mode" "$wall" "$p1" "$p2" "$rss" "$uniq"
        ;;
    esac
}

# -- main loop: one pass per tool -------------------------------------------
echo ""
echo "[bench] start - $(date)"
for pass in "${PASSES[@]}"; do
    set -- $pass
    tool="$1"; mode="$2"

    todo=()
    for c in "${COVERAGES[@]}"; do
        have_run "$c" "$tool" "$mode" || todo+=("$c")
    done
    echo ""
    echo "################ pass $tool $mode - ${#todo[@]}/${#COVERAGES[@]} levels to measure"
    if (( ${#todo[@]} == 0 )); then
        echo "  already complete, skipped"
        continue
    fi

    prev=""
    for c in "${todo[@]}"; do
        want_bytes=$(( CHUNK_BYTES * c ))
        if (( want_bytes >= STOP_BYTES )); then
            echo "  [STOP] ${c}x would exceed $(numfmt --to=iec "$STOP_BYTES")"
            break
        fi
        avail=$(free_bytes "$DATA")
        if (( want_bytes * 3 / 2 > avail )); then
            echo "  [STOP] ${c}x needs ~$(numfmt --to=iec "$want_bytes"), only $(numfmt --to=iec "$avail") free"
            break
        fi
        printf "  %sx  (%s)\n" "$c" "$(date +%H:%M:%S)"
        if ! build_level "$c" "$prev"; then
            echo "  [STOP] failed to build ${c}x"
            rm -f "$DATA/cov_${c}x.fastq.gz"
            break
        fi
        prev="$c"
        echo "    size $(numfmt --to=iec "$(stat -c%s "$DATA/cov_${c}x.fastq.gz")")"
        run_one "$c" "$tool" "$mode"
        rm -rf "$WORK"/* 2>/dev/null
    done

    # The last level of the pass has no successor to consume it.
    [[ -n "$prev" ]] && rm -f "$DATA/cov_${prev}x.fastq.gz"
done

echo ""
echo "[bench] done - $(date)"
echo "[bench] $CSV   ($(( $(wc -l < "$CSV") - 1 )) rows)"
[[ -f "$SCSV" ]] && echo "[bench] $SCSV  ($(( $(wc -l < "$SCSV") - 1 )) rows)  + histograms in $STATS/"
