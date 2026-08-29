#!/usr/bin/env bash
# bench_big_data.sh - the two largest read sets, whole.
#
# gallus and human3 are counted over their complete fof. These are the runs
# that decide whether a tool is usable at all on real sequencing volumes, and
# they dominate the total cost of the benchmark suite - run this last.
#
# Tools: tuna, KMC, FastK.
#
# Output: $ROOT/big_data.csv   (key = number of input files)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${TIMEOUT_S:=21600}"          # 6 h; these runs legitimately take hours
source "$HERE/bench_common.sh"

# name : fof : tuna m : kmc format flag
DATASETS=(
    "gallus:$DATA_ROOT/dataset_reads_gallus/fof.list:21:-fq"
    "human3:$DATA_ROOT/dataset_reads_human3/fof.list:21:-fq"
)

bench_init big_data "tuna kmc fastk"

for spec in "${DATASETS[@]}"; do
    IFS=: read -r ds fof m fmt <<< "$spec"
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    n=$(wc -l < "$fof"); M="$m"
    echo ""
    echo "== $ds  (m=$M, $fmt, $n files) =="
    fkf="$AUX/fastkfof_${ds}.list"; fastk_fof "$fof" "$n" "$fkf"
    run_tuna  "$ds" "$n" "" "@$fof"
    run_kmc   "$ds" "$n" "" "@$fof" "$fmt"
    run_fastk "$ds" "$n" "" "$fkf"
    rm -f "$fkf"
done

bench_done
