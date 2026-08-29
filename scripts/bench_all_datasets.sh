#!/usr/bin/env bash
# bench_all_datasets.sh - per-file measurements across the seven datasets.
#
# For each dataset, the first N files of its fof are counted one at a time.
# This is the "typical single input" view: it says what one genome, one
# metagenome sample or one read set costs, rather than how the tools scale.
#
# Tools: tuna and KMC. FastK is not run here - it has no frozen baseline for
# this experiment and fails on a large share of these inputs.
#
# Output: $ROOT/all_datasets.csv   (key = file index, label = file name)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/bench_common.sh"

# name : fof : max files : tuna m : kmc format flag
DATASETS=(
    "ecoli:$DATA_ROOT/dataset_genome_ecoli/fof.list:100:21:-fm"
    "human:$DATA_ROOT/dataset_genome_human/fof.list:10:21:-fm"
    "salmonella:$DATA_ROOT/dataset_pangenome_salmonella/fof.list:100:21:-fm"
    "gut:$DATA_ROOT/dataset_metagenome_gut/fof.list:100:21:-fm"
    "tara:$DATA_ROOT/dataset_metagenome_tara/fof.list:10:21:-fq"
)

bench_init all_datasets "tuna kmc"

for spec in "${DATASETS[@]}"; do
    IFS=: read -r ds fof maxn m fmt <<< "$spec"
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    M="$m"
    (( MAX_FILES > 0 && MAX_FILES < maxn )) && maxn=$MAX_FILES                                  # per-dataset minimizer length
    echo ""
    echo "== $ds  (m=$M, $fmt, first $maxn files) =="
    i=0
    while IFS= read -r f && (( i < maxn )); do
        [[ -z "$f" || ! -f "$f" ]] && continue
        i=$((i+1))
        echo "  [$i] $(basename "$f")"
        run_tuna "$ds" "$i" "$(basename "$f")" "$f"
        run_kmc  "$ds" "$i" "$(basename "$f")" "$f" "$fmt"
    done < "$fof"
done

bench_done
