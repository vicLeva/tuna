#!/usr/bin/env bash
# bench_thread_sweep.sh - thread sweep, T = 1 .. 32.
#
# One fixed input per dataset, counted at increasing thread counts, to show how
# each tool actually uses the cores it is given.
#
# Both input regimes are swept here, because they stress a counter in opposite
# ways and a tool can scale well on one and badly on the other:
#
#   assemblies  error-free and non-redundant, so nearly every k-mer is distinct
#               and seen once or twice. Passed as a fof, so parallelism comes
#               from spreading FILES across workers.
#   reads       roughly half the k-mers are sequencing errors seen exactly once,
#               which inflates the table with junk and shifts the cost into
#               counting. Passed as ONE compressed file, so parallelism has to
#               come from inside a single stream - a different code path.
#
# The dataset column tells the two apart (`ecoli` vs `ecoli_reads`).
#
# KFC took part in this sweep in earlier rounds and has been dropped.
#
# Tools: tuna, KMC, FastK.
#
# Output: $ROOT/thread_sweep.csv   (key = thread count)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/bench_common.sh"

THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}

# name : input (fof or single file) : tuna m : kmc format flag
DATASETS=(
    "ecoli:$DATA_ROOT/dataset_genome_ecoli/fof_1000.list:21:-fm"
    "human:$DATA_ROOT/dataset_genome_human/fof_10.list:21:-fm"
    "ecoli_reads:$DATA_ROOT/data_sequencing/SRR2584863_1.fastq.gz:21:-fq"
    "human_reads:$DATA_ROOT/data_sequencing/SRR622461_1.fastq.gz:21:-fq"
)

bench_init thread_sweep "tuna kmc fastk"

for spec in "${DATASETS[@]}"; do
    IFS=: read -r ds input m fmt <<< "$spec"
    [[ -e "$input" ]] || { echo "  [skip] $ds: no input at $input"; continue; }
    M="$m"
    # A fof is handed to tuna/KMC with an @ prefix; a plain file is not.
    if [[ "$input" == *.list ]]; then
        tk_in="@$input"
        fkf="$AUX/fastkfof_${ds}.list"; fastk_fof "$input" "$(wc -l < "$input")" "$fkf"
    else
        tk_in="$input"
        printf '%s\n' "$input" > "$AUX/one_${ds}.list"
        fkf="$AUX/fastkfof_${ds}.list"; fastk_fof "$AUX/one_${ds}.list" 1 "$fkf"
    fi
    echo ""
    echo "== $ds  (m=$M, $fmt) =="
    for t in $THREAD_LIST; do
        echo "  T=$t"
        THREADS="$t"
        run_tuna  "$ds" "$t" "" "$tk_in"
        run_kmc   "$ds" "$t" "" "$tk_in" "$fmt"
        run_fastk "$ds" "$t" "" "$fkf"
    done
    rm -f "$fkf" "$AUX/one_${ds}.list"
done

bench_done
