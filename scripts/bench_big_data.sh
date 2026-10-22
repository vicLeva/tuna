#!/usr/bin/env bash
# bench_big_data.sh - the two largest read sets, whole.
#
# gallus and human3 are counted over their complete fof. These are the runs
# that decide whether a tool is usable at all on real sequencing volumes, and
# they dominate the total cost of the benchmark suite - run this last.
#
# One invocation measures ONE configuration: the current K and RAM_GB. The
# paper's table reports four of them, two k values against two memory budgets,
# so run it once per configuration with its own ROOT, or let
# run_paper_suite.sh do that. A -DFIXED_K build only accepts its own k, so
# export TUNA_K27 and TUNA_K31 when covering both.
#
# Tools: tuna, KMC, FastK.
#
# Output: $ROOT/big_data.csv   (key = number of input files)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${TIMEOUT_S:=21600}"          # 6 h; these runs legitimately take hours
source "$HERE/bench_common.sh"

# name : fof : kmc format flag
# m comes from M (or FORCE_M), not from the spec.
DATASETS=(
    "gallus:$DATA_ROOT/dataset_reads_gallus/fof.list:-fq"
    "human3:$DATA_ROOT/dataset_reads_human3/fof.list:-fq"
)

bench_init big_data "tuna kmc fastk"

for spec in "${DATASETS[@]}"; do
    IFS=: read -r ds fof fmt <<< "$spec"
    dataset_enabled "$ds" || continue
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    n=$(wc -l < "$fof")
    echo ""
    echo "== $ds  (k=$K, m=${FORCE_M:-auto}, ram=${RAM_GB}GB, $fmt, $n files) =="
    fkf="$AUX/fastkfof_${ds}.list"; fastk_fof "$fof" "$n" "$fkf"
    run_tuna  "$ds" "$n" "" "@$fof"
    run_kmc   "$ds" "$n" "" "@$fof" "$fmt"
    run_fastk "$ds" "$n" "" "$fkf"
    rm -f "$fkf"
done

bench_done
