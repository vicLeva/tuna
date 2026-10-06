#!/usr/bin/env bash
# bench_msweep.sh - minimizer-length (-m) sweep, tuna only.
#
# Picks a default m by measuring, across data types and collection sizes rather
# than on one dataset. m sets how many kmers a superkmer spans, on average
# s = (k-m+2)/2, so it trades record count against record length:
#
#   small m -> long superkmers, fewer records, less header and overlap overhead
#              per base, but fewer distinct minimizers to partition by
#   large m -> short superkmers, more records, more bytes per base, finer
#              partitioning
#
# At k=31 the phase-1 byte cost per input base runs from about 1.0 at m=9 to
# 2.5 at m=25, a factor of two and a half, so m moves memory as much as time.
# That also means m feeds the in-memory/disk decision: the sizes below are
# chosen so every m stays in memory at RAM_GB=256, keeping the comparison
# within one pipeline. The `spilled` column exists to verify that held.
#
# A build made with -DFIXED_K=31 alone carries m in {9,11,...,29}, so the whole
# sweep runs on one binary. Values outside that set need their own build.
#
# One tuna binary per invocation, as everywhere here: to compare branches, run
# once per build with a different TUNA and ROOT and diff the CSVs.
#
# Resumable: a (dataset, n_files, m) triple already in the CSV is skipped, so an
# interrupted sweep can be restarted and will pick up where it stopped.
#
# Output: $ROOT/msweep.csv   (key = dataset + n_files + m)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${ROOT:=/WORKS/vlevallois/expes_tuna/msweep}"
source "$HERE/bench_common.sh"

# name : fof : how many of its files to use
# Chosen to span assembled vs sequenced data, small vs large collections, and
# low vs high redundancy, while staying in memory for every m.
: "${DATASETS:=\
ecoli_200:$DATA_ROOT/dataset_genome_ecoli/fof.list:200 \
ecoli_3500:$DATA_ROOT/dataset_genome_ecoli/fof.list:3500 \
gut_100:$DATA_ROOT/dataset_genome_gut/fof.list:100 \
human_10:$DATA_ROOT/dataset_genome_human/fof.list:10 \
gallus_12:$DATA_ROOT/dataset_reads_gallus/fof.list:12 \
tara_10:$DATA_ROOT/dataset_reads_tara/fof.list:10}"

: "${MS:=9 11 13 15 17 19 21 23 25}"

CSV="$ROOT/msweep.csv"
AUX="$ROOT/aux"
WORK="$ROOT/work"
mkdir -p "$AUX" "$WORK"

err=0
[[ -z "$TUNA" ]]                 && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
(( err )) && exit 1

[[ -f "$CSV" ]] || echo "dataset,n_files,m,n_parts,wall_s,phase1_s,phase2_s,rss_mb,superkmers,unique_kmers,total_kmers,spilled,status" > "$CSV"

echo "[bench] experiment : msweep"
echo "[bench] k=$K threads=$THREADS ram=${RAM_GB}GB timeout=${TIMEOUT_S}s"
echo "[bench] m values   : $MS"
echo "[bench] tuna       : $TUNA ($(stat -c %y "$TUNA" 2>/dev/null | cut -d. -f1))"
echo "[bench] results    : $CSV"

have_run() {
    awk -F, -v d="$1" -v n="$2" -v m="$3" \
        'NR>1 && $1==d && $2==n && $3==m {f=1} END{exit !f}' "$CSV"
}

for spec in $DATASETS; do
    IFS=: read -r ds fof nf <<< "$spec"
    dataset_enabled "$ds" || continue
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }

    avail=$(wc -l < "$fof")
    (( nf > avail )) && { echo "  [skip] $ds: wants $nf files, fof has $avail"; continue; }
    sub="$AUX/subfof_${ds}.list"; head -n "$nf" "$fof" > "$sub"

    echo ""
    echo "== $ds  ($nf files) =="
    for mm in $MS; do
        if have_run "$ds" "$nf" "$mm"; then
            echo "    m=$mm already measured"
            continue
        fi
        tag="${ds}_m${mm}"; tf="$AUX/$tag.time"; se="$AUX/$tag.stderr"
        rm -rf "$WORK/t"; mkdir -p "$WORK/t"
        /usr/bin/time -v -o "$tf" timeout "$TIMEOUT_S" \
            "$TUNA" -k "$K" -m "$mm" -t "$THREADS" -ram "$RAM_GB" -hp \
            -w "$WORK/t/" "@$sub" "$WORK/out.kff" >/dev/null 2>"$se"
        rc=$?; st=$(status_of "$rc")
        rm -f "$WORK/out.kff"; rm -rf "$WORK/t"

        if [[ $st != ok ]]; then
            echo "$ds,$nf,$mm,,,,,,,,,,$st" >> "$CSV"
            echo "    [$st] m=$mm"
            continue
        fi
        echo "$ds,$nf,$mm,$(se_val n_parts "$se"),$(wall_of "$tf"),$(se_val phase1 "$se"),$(se_val phase2 "$se"),$(rss_of "$tf"),$(se_val superkmers "$se"),$(se_val unique_kmers "$se"),$(se_val total_kmers "$se"),$(se_val spilled "$se"),ok" >> "$CSV"
        printf "    m=%-3s parts=%-7s wall=%9ss  p1=%9ss  p2=%9ss  RSS=%8sMB  sk=%s\n" \
            "$mm" "$(se_val n_parts "$se")" "$(wall_of "$tf")" "$(se_val phase1 "$se")" \
            "$(se_val phase2 "$se")" "$(rss_of "$tf")" "$(se_val superkmers "$se")"
    done
    rm -f "$sub"
done

echo ""
echo "[bench] done - $(date)"
echo "[bench] $CSV   ($(( $(wc -l < "$CSV") - 1 )) rows)"
