#!/usr/bin/env bash
# bench_nsweep.sh - partition-count (-n) sweep, tuna only.
#
# Runs one tuna binary over n = 2^5 .. 2^21 (plus the auto-tuned value), on
# whichever datasets are configured, to find how end-to-end time depends on
# the partition count and where that optimum sits relative to auto-tuning.
#
# -n is rounded UP to a power of two by tuna (partition_fn masks with n-1), so
# these are the entire reachable search space at each step; nothing finer
# exists between two consecutive powers of two.
#
# One tuna binary per invocation, like every other script here: to compare
# branches (say kache-hash vs unordered_dense), run this once per build with a
# different TUNA and ROOT each time and diff the CSVs afterwards.
#
# ---------------------------------------------------------------------------
# Two passes, two CSVs
# ---------------------------------------------------------------------------
# pass 1  timing        every n in $NS      -> $ROOT/nsweep.csv
# pass 2  table stats   every n in $DBG_NS  -> $ROOT/nsweep_parts.csv
#
# $DBG_NS defaults to $NS, so every point is measured both ways and each n
# costs TWO full runs of the workload. Set DBG_NS to a subset to sample the
# statistics more coarsely (they move smoothly between points), or to the
# empty string to skip pass 2 entirely.
#
# Pass 2 runs with -dbg, which makes tuna emit two per-run CSVs describing
# every partition's table (debug_table_stats.csv, debug_resize_events.csv).
# That bookkeeping perturbs wall time, so pass 2's own timings are NOT
# comparable to pass 1's - it exists to look inside the tables, not to race
# them, which is why it is a separate run rather than -dbg on the timed one.
#
# Table layout differs by branch (kache-hash quotients the key out of the
# bucket address; unordered_dense and absl::flat_hash_map store it in full),
# so "table_cap" here means whatever that branch calls capacity/bucket_count -
# converting it to bytes/partition is left to analysis, per branch, not to
# this script. On a branch whose table does not track resizes (e.g.
# unordered_dense's resize_log() is a stub), n_resizes is always 0: that is a
# property of the table, not a measurement failure.
#
# Output: $ROOT/nsweep.csv         (key = n)              - pass 1
#         $ROOT/nsweep_parts.csv   (key = n)              - pass 2
#         $ROOT/dbg/<dataset>_n<n>/                        - pass 2 raw CSVs

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${ROOT:=/WORKS/vlevallois/expes_tuna/nsweep}"
source "$HERE/bench_common.sh"

# name : fof : tuna m
: "${DATASETS:=gallus:$DATA_ROOT/dataset_reads_gallus/fof.list:21 human3:$DATA_ROOT/dataset_reads_human3/fof.list:21}"
: "${NS:=auto 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152}"
: "${DBG_NS:=$NS}"          # every n by default; narrow it to sample fewer points

CSV="$ROOT/nsweep.csv"
PCSV="$ROOT/nsweep_parts.csv"
DBG="$ROOT/dbg"
AUX="$ROOT/aux/nsweep"
WORK="$ROOT/work"
mkdir -p "$DBG" "$AUX" "$WORK"

err=0
[[ -z "$TUNA" ]]              && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
(( err )) && exit 1

[[ -f "$CSV" ]]  || echo "dataset,n,n_parts,wall_s,phase1_s,phase2_s,rss_mb,unique_kmers,total_kmers,status" > "$CSV"
[[ -f "$PCSV" ]] || echo "dataset,n,n_parts,parts_nonempty,uniq_mean,uniq_median,uniq_min,uniq_max,inserted_mean,init_sz_mean,init_cap_mean,cap_mean,load_mean,parts_grown,n_resizes_total,parts_resized" > "$PCSV"

echo "[bench] experiment : nsweep"
echo "[bench] k=$K m=$M threads=$THREADS ram=${RAM_GB}GB timeout=${TIMEOUT_S}s"
echo "[bench] tuna       : $TUNA ($(stat -c %y "$TUNA" 2>/dev/null | cut -d. -f1))"
echo "[bench] timings    : $CSV"
echo "[bench] table stats: $PCSV"

have_run()  { awk -F, -v d="$1" -v n="$2" 'NR>1 && $1==d && $2==n {f=1} END{exit !f}' "$CSV"; }
have_stat() { awk -F, -v d="$1" -v n="$2" 'NR>1 && $1==d && $2==n {f=1} END{exit !f}' "$PCSV"; }
in_dbg_ns() { case " $DBG_NS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

# args: <dataset> <fof> <m> <n> <dbg flag: "" or "-dbg"> <out dir: "" or path>
# Fills globals: WALL P1 P2 NP UNIQ TOT ST
run_once() {
    local ds="$1" fof="$2" mm="$3" n="$4" dbg="$5" out="$6" nopt=""
    [[ "$n" != auto ]] && nopt="-n $n"
    local tag="${ds}_n${n}${dbg:+_dbg}"; local tf="$AUX/$tag.time" se="$AUX/$tag.stderr"
    rm -rf "$WORK/t"; mkdir -p "$WORK/t"
    /usr/bin/time -v -o "$tf" timeout "$TIMEOUT_S" \
        "$TUNA" -k "$K" -m "$mm" -t "$THREADS" -ram "$RAM_GB" $nopt -hp $dbg \
        -w "$WORK/t/" "@$fof" "$WORK/out.kff" >/dev/null 2>"$se"
    local rc=$?; ST=$(status_of "$rc")
    if [[ -n "$out" ]]; then
        mkdir -p "$out"
        cp -f "$WORK/t/debug_table_stats.csv"   "$out/" 2>/dev/null
        cp -f "$WORK/t/debug_resize_events.csv" "$out/" 2>/dev/null
        cp -f "$se"                             "$out/stderr.txt" 2>/dev/null
    fi
    rm -f "$WORK/out.kff"; rm -rf "$WORK/t"
    WALL=$(wall_of "$tf"); P1=$(se_val phase1 "$se"); P2=$(se_val phase2 "$se")
    NP=$(se_val n_parts "$se"); UNIQ=$(se_val unique_kmers "$se"); TOT=$(se_val total_kmers "$se")
    RSS=$(rss_of "$tf")
}

for spec in $DATASETS; do
    IFS=: read -r ds fof mm <<< "$spec"
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    mm="${FORCE_M:-$mm}"
    echo ""
    echo "== $ds  (m=$mm) =="

    for n in $NS; do
        # ---- pass 1: timing ----
        if have_run "$ds" "$n"; then
            echo "    [time]  n=$n already measured"
        else
            run_once "$ds" "$fof" "$mm" "$n" "" ""
            if [[ $ST != ok ]]; then
                echo "$ds,$n,,,,,,,,$ST" >> "$CSV"
                echo "    [time]  [$ST] n=$n"
            else
                echo "$ds,$n,$NP,$WALL,$P1,$P2,$RSS,$UNIQ,$TOT,ok" >> "$CSV"
                printf "    [time]  n=%-9s parts=%-9s wall=%9ss  p1=%9ss  p2=%9ss  RSS=%8sMB\n" \
                    "$n" "$NP" "$WALL" "$P1" "$P2" "$RSS"
            fi
        fi

        # ---- pass 2: table stats (every n in DBG_NS) ----
        in_dbg_ns "$n" || continue
        if have_stat "$ds" "$n"; then
            echo "    [stats] n=$n already measured"
            continue
        fi
        out="$DBG/${ds}_n${n}"
        run_once "$ds" "$fof" "$mm" "$n" "-dbg" "$out"
        if [[ $ST != ok || ! -s "$out/debug_table_stats.csv" ]]; then
            echo "$ds,$n,,,,,,,,,,,,,," >> "$PCSV"
            echo "    [stats] [$ST] n=$n"
            continue
        fi
        # partition_id,init_sz,init_cap,table_cap,n_inserted,n_unique,load_factor,n_resizes,resize_s
        # median via sort: asort() is a gawk extension and the node may run mawk
        med=$(awk -F, 'NR>1 && $6>0 {print $6}' "$out/debug_table_stats.csv" | sort -n | \
              awk '{a[NR]=$1} END{ if(NR==0) print 0;
                                   else if(NR%2) print a[(NR+1)/2];
                                   else print (a[NR/2]+a[NR/2+1])/2 }')
        awk -F, -v ds="$ds" -v n="$n" -v np="$NP" -v med="$med" '
            NR>1 {
                parts++; if($6>0) ne++
                su+=$6; si+=$5; sinit+=$2; sicap+=$3; scap+=$4; sload+=$7; sres+=$8
                # table_cap > init_cap means the table outgrew the reserve
                # and rehashed - the only way to see that on a table keeping
                # no resize log, which is every table except kache-hash.
                if($4>$3) grew++
                if($8>0) pres++
                if(mn=="" || $6<mn) mn=$6
                if($6>mx) mx=$6
            }
            END{
                if(ne==0){ print ds","n","np",0,,,,,,,,,,,"; exit }
                printf "%s,%s,%s,%d,%.1f,%.1f,%d,%d,%.1f,%.1f,%.1f,%.1f,%.4f,%d,%d,%d\n",
                       ds,n,np,ne,su/ne,med,mn,mx,si/ne,sinit/ne,sicap/ne,scap/ne,sload/ne,grew+0,sres,pres+0
            }' "$out/debug_table_stats.csv" >> "$PCSV"
        printf "    [stats] n=%-9s %s\n" "$n" \
            "$(tail -1 "$PCSV" | awk -F, '{printf "parts=%s uniq mean/max=%.0f/%s cap_mean=%.0f load=%.3f grew=%s",$4,$5,$8,$12,$13,$14}')"
    done
done

echo ""
echo "[bench] done - $(date)"
echo "[bench] $CSV   ($(( $(wc -l < "$CSV") - 1 )) rows)"
echo "[bench] $PCSV  ($(( $(wc -l < "$PCSV") - 1 )) rows)"
echo "[bench] per-run debug CSVs under $DBG/"
