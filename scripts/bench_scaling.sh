#!/usr/bin/env bash
# bench_scaling.sh - how the tools scale with the number of input files.
#
# Counts the first n files of a fof, for increasing n. Two collections:
# assembled E. coli genomes (many small files) and human assemblies (few large
# ones), which stress the tools in opposite ways.
#
# Tools: tuna, KMC, FastK.
#
# FastK is given a shorter timeout than the others: it is known to hang rather
# than fail on certain file-count/thread-count combinations, and there is no
# point spending six hours confirming it. Hangs land in the CSV as `timeout`,
# crashes as `fail`.
#
# Output: $ROOT/scaling.csv   (key = number of input files)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${FASTK_TIMEOUT_S:=900}"
source "$HERE/bench_common.sh"

ECOLI_FOF="$DATA_ROOT/dataset_genome_ecoli/fof.list"
HUMAN_FOF="$DATA_ROOT/dataset_genome_human/fof.list"
ECOLI_NS=(1 2 3 5 10 20 50 100 200 500 1000 1500 2000 2500 3000 3500)
HUMAN_NS=(1 2 3 4 5 6 7 8 9 10 15 20 25 30 60)

bench_init scaling "tuna kmc fastk"

sweep() {
    local ds="$1" fof="$2" m="$3"; shift 3
    local ns=("$@") total n sub fkf
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; return; }
    total=$(wc -l < "$fof")
    M="$m"
    echo ""
    echo "== $ds  (m=$M, $total files available) =="
    for n in "${ns[@]}"; do
        (( n > total )) && { echo "  n=$n > $total, skipped"; continue; }
        echo "  n=$n"
        sub="$AUX/subfof_${ds}_${n}.list"; head -n "$n" "$fof" > "$sub"
        fkf="$AUX/fastkfof_${ds}_${n}.list"; fastk_fof "$fof" "$n" "$fkf"
        run_tuna  "$ds" "$n" "" "@$sub"
        run_kmc   "$ds" "$n" "" "@$sub" -fm
        run_fastk "$ds" "$n" "" "$fkf"
        rm -f "$sub" "$fkf"
    done
}

sweep ecoli "$ECOLI_FOF" 21 "${ECOLI_NS[@]}"
sweep human "$HUMAN_FOF" 21 "${HUMAN_NS[@]}"

bench_done
