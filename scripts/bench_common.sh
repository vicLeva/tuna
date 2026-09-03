#!/usr/bin/env bash
# bench_common.sh - shared machinery for the tuna paper benchmarks.
#
# Sourced by every bench_*.sh in this directory. Not runnable on its own.
#
# ---------------------------------------------------------------------------
# What is measured
# ---------------------------------------------------------------------------
# Two output modes per tool, so the cost of writing ASCII can be isolated:
#
#   bin    the tool's own binary output   (tuna KFF, KMC .kmc_pre/.kmc_suf,
#                                          FastK .ktab)
#   ascii  a plain-text k-mer/count table (TSV)
#
# The three tools reach ASCII differently, and that shapes how they are timed:
#
#   tuna   writes ASCII directly, in a single pass. It therefore needs TWO
#          separate runs - one producing KFF, one producing TSV - and the
#          ASCII cost is (ascii_wall - bin_wall).
#   KMC    writes a binary DB, then kmc_dump converts it. ONE kmc run plus a
#          separately timed kmc_dump: ascii_wall = kmc_wall + dump_wall.
#   FastK  writes a .ktab table, then Tabex converts it. Same shape as KMC:
#          ascii_wall = fastk_wall + tabex_wall.
#
# So for KMC and FastK the `dump_s` column holds the conversion step on its
# own; for tuna it is empty and the equivalent quantity is ascii - bin.
#
# ---------------------------------------------------------------------------
# CSV schema (one file per experiment)
# ---------------------------------------------------------------------------
#   dataset,key,label,tool,mode,wall_s,rss_mb,phase1_s,phase2_s,dump_s,
#   unique_kmers,total_kmers,status
#
#   key    the experiment's independent variable: file index (all_datasets),
#          number of input files (scaling, big_data) or thread count
#          (thread_sweep).
#   label  free text for that point; the file name in all_datasets, empty
#          elsewhere.
#   phase1_s / phase2_s
#          tuna: its own phase1/phase2. KMC: 1st/2nd stage. FastK: empty,
#          it reports no equivalent breakdown.
#   status ok | fail | timeout. Rows are always written, including failures:
#          "FastK cannot do this" is a result and must survive into the CSV.
#
# ---------------------------------------------------------------------------
# Failure handling
# ---------------------------------------------------------------------------
# Every run is wrapped in `timeout`, because FastK is known to hang outright
# on some file-count/thread-count combinations (it spins in Scan_All_Input and
# never returns). A run killed by the timeout is recorded as `timeout`; any
# other non-zero exit is `fail`. The two are kept distinct: a hang and a crash
# say different things about a tool.
#
# Tabex crashes on a large fraction of FastK tables. When that happens the
# `bin` row still stands and only the `ascii` row is marked failed.

set -uo pipefail
export LC_ALL=C LANG=C

# ---------------------------------------------------------------------------
# Configuration (override from the environment)
# ---------------------------------------------------------------------------
: "${TUNA:=}"                # required
: "${KMC:=}"                 # required
: "${KMC_DUMP:=}"            # required
: "${FASTK:=}"               # required only by experiments that use FastK
: "${TABEX:=}"
: "${KMC_TOOLS:=}"           # optional: enables per-input k-mer statistics

: "${ROOT:=/WORKS/vlevallois/expes_tuna/expes_paper}"
# Root of the input collections. Overridable so the scripts can be exercised
# on a small fake tree without editing them.
: "${DATA_ROOT:=/WORKS/vlevallois/data}"

: "${K:=31}"
: "${M:=21}"                 # tuna minimizer length
: "${THREADS:=8}"
: "${RAM_GB:=256}"

: "${TIMEOUT_S:=21600}"      # 6 h per run, all tools

# --- overrides, mainly for side experiments -------------------------------
# FORCE_M    ignore the per-dataset minimizer length and use this one
#            everywhere (for choosing m by comparing whole experiments).
# ONLY_TOOLS space-separated subset of "tuna kmc fastk" to actually run.
# MAX_FILES  cap on how many files a dataset contributes, where the
#            experiment iterates over files.
: "${FORCE_M:=}"
: "${ONLY_TOOLS:=}"
: "${MAX_FILES:=0}"
# ONLY_MODES  subset of "bin ascii" to measure. Comparing two counting-table
#             implementations only needs `bin`: the ascii run repeats the same
#             counting work and adds TSV formatting, which such a change does
#             not touch, so it doubles the cost for no signal.
: "${ONLY_MODES:=}"
# KEEP_OUTPUT  directory in which to preserve the BINARY output of each run
#              instead of deleting it (KMC .kmc_pre/.kmc_suf, tuna .kff).
#              ASCII output is never kept: it is ~40 bytes per distinct k-mer
#              against ~5 for binary, which runs to hundreds of GB on the big
#              read sets. Everything a k-mer histogram can tell you is already
#              in the binary DB - see scripts/kmer_stats.sh.
: "${KEEP_OUTPUT:=}"
: "${FASTK_TIMEOUT_S:=$TIMEOUT_S}"

WORK="$ROOT/work"; AUX="$ROOT/aux"; LINKS="$ROOT/fastk_links"

# ---------------------------------------------------------------------------
# Setup / validation
# ---------------------------------------------------------------------------
# args: <experiment name> <"tool tool ..."> - the tools this experiment needs
bench_init() {
    EXPERIMENT="$1"; TOOLS_USED="$2"
    CSV="$ROOT/${EXPERIMENT}.csv"
    AUX="$ROOT/aux/$EXPERIMENT"

    # Only demand the binaries that will actually be used.
    [[ -n "$ONLY_TOOLS" ]] && TOOLS_USED="$ONLY_TOOLS"

    local err=0 v
    if [[ "$TOOLS_USED" == *tuna* ]]; then
        [[ -z "$TUNA" ]] && { echo "[error] TUNA is not set"; err=1; }
        [[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
    fi
    if [[ "$TOOLS_USED" == *kmc* ]]; then
        for v in KMC KMC_DUMP; do
            [[ -z "${!v}" ]] && { echo "[error] $v is not set"; err=1; }
            [[ -n "${!v}" && ! -x "${!v}" ]] && { echo "[error] not executable: ${!v}"; err=1; }
        done
    fi
    if [[ "$TOOLS_USED" == *fastk* ]]; then
        for v in FASTK TABEX; do
            [[ -z "${!v}" ]] && { echo "[error] $v is not set (needed by $EXPERIMENT)"; err=1; }
            [[ -n "${!v}" && ! -x "${!v}" ]] && { echo "[error] not executable: ${!v}"; err=1; }
        done
    fi
    (( err )) && exit 1

    mkdir -p "$ROOT" "$WORK" "$AUX" "$LINKS"
    [[ -f "$CSV" ]] || echo "dataset,key,label,tool,mode,wall_s,rss_mb,phase1_s,phase2_s,dump_s,unique_kmers,total_kmers,status" > "$CSV"

    SCSV="$ROOT/${EXPERIMENT}_kmer_stats.csv"
    HIST="$ROOT/hist/$EXPERIMENT"
    if [[ -n "$KMC_TOOLS" ]]; then
        mkdir -p "$HIST"
        [[ -f "$SCSV" ]] || echo "dataset,key,label,distinct_kmers,total_kmers,singletons,pct_singletons,mean_count,median_count,p99_count,max_count,cov_peak,top1pct_mass_pct,gini" > "$SCSV"
        echo "[bench] k-mer stats : $SCSV"
    fi

    echo "[bench] experiment : $EXPERIMENT"
    echo "[bench] tools      : $TOOLS_USED"
    [[ -n "$FORCE_M" ]] && echo "[bench] FORCE_M=$FORCE_M (per-dataset m ignored)"
    [[ "$MAX_FILES" != 0 ]] && echo "[bench] MAX_FILES=$MAX_FILES"
    echo "[bench] k=$K m=$M threads=$THREADS ram=${RAM_GB}GB timeout=${TIMEOUT_S}s (fastk ${FASTK_TIMEOUT_S}s)"
    echo "[bench] tuna       : $TUNA ($(stat -c %y "$TUNA" 2>/dev/null | cut -d. -f1))"
    echo "[bench] output     : $CSV"
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
wall_of() { awk -F: '/Elapsed \(wall clock\)/{t=$0} END{n=split(t,a," ");split(a[n],b,":");
            if(length(b)==3) printf "%.3f",b[1]*3600+b[2]*60+b[3]; else printf "%.3f",b[1]*60+b[2]}' "$1" 2>/dev/null; }
rss_of()  { awk '/Maximum resident/{printf "%.0f",$NF/1024}' "$1" 2>/dev/null; }
se_val()  { grep -m1 "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2;gsub(/^ +/,"",v);sub(/s$/,"",v);print v}'; }
kmc_num() { grep -m1 "$1" "$2" 2>/dev/null | awk -F: '{gsub(/^ +/,"",$2);print $2}'; }

row() { echo "$1" >> "$CSV"; }

# Is this output mode part of the run? (ONLY_MODES empty means both)
mode_enabled() { [[ -z "$ONLY_MODES" || " $ONLY_MODES " == *" $1 "* ]]; }

have_stats() { [[ -f "${SCSV:-}" ]] && awk -F, -v d="$1" -v k="$2" 'NR>1 && $1==d && $2==k {f=1} END{exit !f}' "$SCSV"; }

# ---------------------------------------------------------------------------
# k-mer count statistics
# ---------------------------------------------------------------------------
# Everything is derived from the count histogram rather than from an ASCII
# dump: the histogram is a few hundred kB whatever the input size, while a TSV
# is ~40 bytes per distinct k-mer - hundreds of GB on the big read sets - and a
# pass over it would cost hours for numbers that are already implied.
#
#   singletons     k-mers seen exactly once: mostly sequencing error in reads
#   cov_peak       most common count above 1, i.e. the coverage mode
#   top1pct_mass   share of all occurrences held by the 1% most abundant k-mers
#   gini           inequality of the count distribution, 0 = flat, 1 = extreme
#
# args: <csv field prefix> <kmc db prefix> <histogram path>
# Echoes one CSV line, or nothing if the histogram could not be produced.
kmer_stats_line() {
    local prefix="$1" db="$2" hist="$3"
    [[ -z "$KMC_TOOLS" || ! -x "$KMC_TOOLS" ]] && return 0
    [[ -s "${db}.kmc_suf" ]] || return 0
    "$KMC_TOOLS" transform "$db" histogram "$hist" -ci1 -cx1000000 >/dev/null 2>&1 || return 0
    awk -v pre="$prefix" '
        { c=$1+0; n=$2+0; if(n==0) next
          cnt[++m]=c; num[m]=n; d+=n; t+=c*n; if(c==1) s1=n; if(c>mx) mx=c
          if(c>1 && n>peakn){peakn=n; peak=c} }
        END{
          if(d==0){ exit }
          need50=d*0.5; need99=d*0.99; need1=d*0.01
          cum=0; for(i=1;i<=m;i++){ cum+=num[i]
            if(!med && cum>=need50) med=cnt[i]
            if(!p99 && cum>=need99) p99=cnt[i] }
          cum=0; mass=0; for(i=m;i>=1;i--){ take=num[i]
            if(cum+take>need1) take=need1-cum
            mass+=take*cnt[i]; cum+=take; if(cum>=need1) break }
          cum=0; area=0; for(i=1;i<=m;i++){ area+=num[i]*(2*cum+num[i])*cnt[i]; cum+=num[i] }
          gini=(d>1&&t>0)? area/(d*t)-1 : 0
          printf "%s,%d,%d,%d,%.3f,%.3f,%d,%d,%d,%d,%.2f,%.4f\n",
                 pre,d,t,s1,100*s1/d,t/d,med,p99,mx,peak,100*mass/t,gini
        }' "$hist"
}

# Is this tool part of the run? (ONLY_TOOLS empty means "all of them")
tool_enabled() { [[ -z "$ONLY_TOOLS" || " $ONLY_TOOLS " == *" $1 "* ]]; }

# Resume at measurement granularity: re-running a script only fills the gaps.
have_run() {
    awk -F, -v d="$1" -v k="$2" -v t="$3" -v m="$4" \
        'NR>1 && $1==d && $2==k && $4==t && $5==m {f=1} END{exit !f}' "$CSV"
}

# rc 124 is what `timeout` returns when it had to kill the process.
status_of() { [[ "$1" -eq 0 ]] && echo ok || { [[ "$1" -eq 124 ]] && echo timeout || echo fail; }; }

# ---------------------------------------------------------------------------
# FastK input links
# ---------------------------------------------------------------------------
# FastK dispatches on file extension and does not recognise .fna, so inputs are
# exposed through a symlink farm with names it accepts. Cheap, and it leaves
# the real data untouched.
fastk_link() {
    local f="$1" base; base=$(basename "$f")
    case "$base" in
        *.fna.gz) echo "$LINKS/${base%.fna.gz}.fa.gz" ;;
        *.fna)    echo "$LINKS/${base%.fna}.fa"       ;;
        *)        echo "$LINKS/$base"                 ;;
    esac
}

# args: <fof> <n> <out fof>  - first n entries of fof, as FastK-safe links
fastk_fof() {
    local fof="$1" n="$2" out="$3" f lnk
    : > "$out"
    while IFS= read -r f; do
        [[ -z "$f" ]] && continue
        lnk=$(fastk_link "$f"); ln -sf "$f" "$lnk"; echo "$lnk" >> "$out"
    done < <(head -n "$n" "$fof")
}

# ---------------------------------------------------------------------------
# Tool runners
# ---------------------------------------------------------------------------
# Each writes one bin row and one ascii row, and prints one line per mode.
#
# Common args: <dataset> <key> <label> <input>
#   <input> is what the tool should read: "@/path/fof.list" for a file of
#   files, or a plain path for a single file.

run_tuna() {
    tool_enabled tuna || return 0
    local ds="$1" key="$2" label="$3" input="$4" mode out tag tf se rc st
    local mm="${FORCE_M:-$M}"
    for mode in bin ascii; do
        mode_enabled "$mode" || continue
        have_run "$ds" "$key" tuna "$mode" && { echo "    [tuna $mode] already measured"; continue; }
        tag="${ds}_${key}_tuna_${mode}"; tf="$AUX/$tag.time"; se="$AUX/$tag.stderr"
        [[ $mode == bin ]] && out="$WORK/out.kff" || out="$WORK/out.tsv"
        rm -rf "$WORK/t"; mkdir -p "$WORK/t"
        /usr/bin/time -v -o "$tf" timeout "$TIMEOUT_S" \
            "$TUNA" -k "$K" -m "$mm" -t "$THREADS" -ram "$RAM_GB" -hp \
            -w "$WORK/t/" "$input" "$out" >/dev/null 2>"$se"
        rc=$?; st=$(status_of "$rc")
        if [[ -n "$KEEP_OUTPUT" && $mode == bin && $st == ok ]]; then
            mkdir -p "$KEEP_OUTPUT"
            mv -f "$out" "$KEEP_OUTPUT/${ds}_${key}_tuna.kff" 2>/dev/null
        fi
        rm -f "$out"; rm -rf "$WORK/t"
        if [[ $st != ok ]]; then
            row "$ds,$key,$label,tuna,$mode,,,,,,,,$st"; echo "    [$st] tuna $mode"; continue
        fi
        row "$ds,$key,$label,tuna,$mode,$(wall_of "$tf"),$(rss_of "$tf"),$(se_val phase1 "$se"),$(se_val phase2 "$se"),,$(se_val unique_kmers "$se"),$(se_val total_kmers "$se"),ok"
        printf "    [tuna  %-5s] wall=%9ss  p1=%9ss  p2=%9ss  RSS=%8sMB\n" \
            "$mode" "$(wall_of "$tf")" "$(se_val phase1 "$se")" "$(se_val phase2 "$se")" "$(rss_of "$tf")"
    done
}

# extra arg: <kmc format flag>  (-fm multi-FASTA, -fq FASTQ)
run_kmc() {
    tool_enabled kmc || return 0
    local ds="$1" key="$2" label="$3" input="$4" fmt="$5"
    local tag tf se dt rc st wall rss p1 p2 uniq tot dump drc dst
    local need_bin=1 need_ascii=1
    have_run "$ds" "$key" kmc bin   && need_bin=0
    have_run "$ds" "$key" kmc ascii && need_ascii=0
    mode_enabled ascii || need_ascii=0
    mode_enabled bin   || need_bin=0
    (( need_bin == 0 && need_ascii == 0 )) && { echo "    [kmc] already measured"; return; }

    tag="${ds}_${key}_kmc"; tf="$AUX/$tag.time"; se="$AUX/$tag.stderr"; dt="$AUX/$tag.dumptime"
    rm -rf "$WORK/kmctmp"; mkdir -p "$WORK/kmctmp"
    /usr/bin/time -v -o "$tf" timeout "$TIMEOUT_S" \
        "$KMC" -k"$K" -ci1 -cs4294967295 "$fmt" -m"$RAM_GB" -hp -t"$THREADS" \
        "$input" "$WORK/kmcdb" "$WORK/kmctmp" > "$se" 2>&1
    rc=$?; st=$(status_of "$rc"); rm -rf "$WORK/kmctmp"
    if [[ $st != ok ]]; then
        (( need_bin ))   && row "$ds,$key,$label,kmc,bin,,,,,,,,$st"
        (( need_ascii )) && row "$ds,$key,$label,kmc,ascii,,,,,,,,$st"
        echo "    [$st] kmc"; rm -f "$WORK/kmcdb".kmc_*; return
    fi
    wall=$(wall_of "$tf"); rss=$(rss_of "$tf")
    p1=$(se_val "1st stage" "$se"); p2=$(se_val "2nd stage" "$se")
    uniq=$(kmc_num "No. of unique k-mers" "$se"); tot=$(kmc_num "Total no. of k-mers" "$se")
    if (( need_bin )); then
        row "$ds,$key,$label,kmc,bin,$wall,$rss,$p1,$p2,,$uniq,$tot,ok"
        printf "    [kmc   %-5s] wall=%9ss  p1=%9ss  p2=%9ss  RSS=%8sMB\n" bin "$wall" "$p1" "$p2" "$rss"
    fi
    # ASCII = the same binary DB, plus a separately timed kmc_dump.
    if (( need_ascii )); then
        /usr/bin/time -v -o "$dt" timeout "$TIMEOUT_S" \
            "$KMC_DUMP" "$WORK/kmcdb" "$WORK/out.tsv" >/dev/null 2>>"$se"
        drc=$?; dst=$(status_of "$drc"); rm -f "$WORK/out.tsv"
        if [[ $dst != ok ]]; then
            row "$ds,$key,$label,kmc,ascii,,,,,,,,$dst"; echo "    [$dst] kmc_dump"
        else
            dump=$(wall_of "$dt")
            row "$ds,$key,$label,kmc,ascii,$(awk "BEGIN{printf \"%.3f\",$wall+$dump}"),$rss,$p1,$p2,$dump,$uniq,$tot,ok"
            printf "    [kmc   %-5s] wall=%9ss  (kmc=%ss + dump=%ss)\n" ascii \
                "$(awk "BEGIN{printf \"%.3f\",$wall+$dump}")" "$wall" "$dump"
        fi
    fi
    # Statistics come from the database this run just built, before it goes.
    if [[ -n "$KMC_TOOLS" ]] && ! have_stats "$ds" "$key"; then
        local line; line=$(kmer_stats_line "$ds,$key,$label" "$WORK/kmcdb" "$HIST/${ds}_${key}.hist")
        if [[ -n "$line" ]]; then
            echo "$line" >> "$SCSV"
            echo "    [stats] $(echo "$line" | awk -F, '{printf "distinct=%s singletons=%.1f%% peak=%sx top1%%=%s%% gini=%s",$4,$7,$12,$13,$14}')"
        fi
    fi
    if [[ -n "$KEEP_OUTPUT" ]]; then
        mkdir -p "$KEEP_OUTPUT"
        mv -f "$WORK/kmcdb.kmc_pre" "$KEEP_OUTPUT/${ds}_${key}.kmc_pre" 2>/dev/null
        mv -f "$WORK/kmcdb.kmc_suf" "$KEEP_OUTPUT/${ds}_${key}.kmc_suf" 2>/dev/null
        echo "    [keep] $KEEP_OUTPUT/${ds}_${key}.kmc_{pre,suf}  ($(du -sh "$KEEP_OUTPUT" 2>/dev/null | cut -f1) kept so far)"
    fi
    rm -f "$WORK/kmcdb".kmc_*
}

# args: <dataset> <key> <label> <fastk fof of symlinks>
run_fastk() {
    tool_enabled fastk || return 0
    local ds="$1" key="$2" label="$3" fof="$4"
    local tag tf se dt rc st wall rss dump drc dst db tmp
    local need_bin=1 need_ascii=1
    have_run "$ds" "$key" fastk bin   && need_bin=0
    have_run "$ds" "$key" fastk ascii && need_ascii=0
    mode_enabled ascii || need_ascii=0
    mode_enabled bin   || need_bin=0
    (( need_bin == 0 && need_ascii == 0 )) && { echo "    [fastk] already measured"; return; }

    tag="${ds}_${key}_fastk"; tf="$AUX/$tag.time"; se="$AUX/$tag.stderr"; dt="$AUX/$tag.dumptime"
    db="$WORK/fastk_$tag"; tmp="$WORK/fastktmp_$tag"
    rm -rf "$tmp"; mkdir -p "$tmp"
    mapfile -t files < "$fof"
    # -t1 makes FastK write the .ktab table; without it only a histogram is
    # produced and there would be nothing for Tabex to convert.
    /usr/bin/time -v -o "$tf" timeout "$FASTK_TIMEOUT_S" \
        "$FASTK" -k"$K" -t1 -T"$THREADS" -M"$RAM_GB" -N"$db" -P"$tmp" "${files[@]}" \
        >/dev/null 2>"$se"
    rc=$?; st=$(status_of "$rc"); rm -rf "$tmp"
    if [[ $st != ok ]]; then
        (( need_bin ))   && row "$ds,$key,$label,fastk,bin,,,,,,,,$st"
        (( need_ascii )) && row "$ds,$key,$label,fastk,ascii,,,,,,,,$st"
        echo "    [$st] fastk"; rm -f "$db"* ".${db##*/}"* 2>/dev/null; return
    fi
    wall=$(wall_of "$tf"); rss=$(rss_of "$tf")
    if (( need_bin )); then
        row "$ds,$key,$label,fastk,bin,$wall,$rss,,,,,,ok"
        printf "    [fastk %-5s] wall=%9ss  RSS=%8sMB\n" bin "$wall" "$rss"
    fi
    # Tabex crashes on a good share of FastK tables; when it does, only the
    # ascii row fails and the bin measurement above still stands.
    if (( need_ascii )); then
        /usr/bin/time -v -o "$dt" timeout "$FASTK_TIMEOUT_S" \
            "$TABEX" -A "$db" > "$WORK/out.tsv" 2>>"$se"
        drc=$?; dst=$(status_of "$drc"); rm -f "$WORK/out.tsv"
        if [[ $dst != ok ]]; then
            row "$ds,$key,$label,fastk,ascii,,,,,,,,$dst"; echo "    [$dst] Tabex"
        else
            dump=$(wall_of "$dt")
            row "$ds,$key,$label,fastk,ascii,$(awk "BEGIN{printf \"%.3f\",$wall+$dump}"),$rss,,,$dump,,,ok"
            printf "    [fastk %-5s] wall=%9ss  (fastk=%ss + tabex=%ss)\n" ascii \
                "$(awk "BEGIN{printf \"%.3f\",$wall+$dump}")" "$wall" "$dump"
        fi
    fi
    rm -f "$db"* 2>/dev/null
    rm -f "$(dirname "$db")/.$(basename "$db")"* 2>/dev/null
}

bench_done() {
    echo ""
    echo "[bench] done - $(date)"
    echo "[bench] $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
    awk -F, 'NR>1{print "         "$4"-"$5"-"$13}' "$CSV" | sort | uniq -c
}
