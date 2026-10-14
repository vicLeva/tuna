#!/usr/bin/env bash
# bench_all_datasets.sh - per-file and whole-collection measurements.
#
# Two passes over the same files, for each dataset:
#
#   per file     the first N files counted one at a time. The "typical single
#                input" view: what one genome, one metagenome sample or one
#                read set costs.
#   whole        those same N files counted together in a single run. What a
#                user pays to count the collection in one go, which is the
#                regime where peak memory and the partitioning cost matter.
#
# The whole-collection pass takes exactly the files the per-file pass used, so
# the two views describe the same data and can be compared directly.
#
# Tools: tuna, KMC and FastK. FastK fails on a large share of these inputs
# (Tabex segfaults on the big ones), so expect fail rows rather than results
# on human and tara; the bin measurement still stands where only ascii fails.
#
# Output: $ROOT/all_datasets.csv
#   key = file index, label = file name   for the per-file rows
#   key = "all",      label = file count  for the whole-collection rows

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/bench_common.sh"

# name : fof : max files : kmc format flag
# m is not per dataset: it comes from M (or FORCE_M) so the whole suite
# measures one minimizer length.
DATASETS=(
    "ecoli:$DATA_ROOT/dataset_genome_ecoli/fof.list:100:-fm"
    "human:$DATA_ROOT/dataset_genome_human/fof.list:10:-fm"
    "salmonella:$DATA_ROOT/dataset_pangenome_salmonella/fof.list:100:-fm"
    "gut:$DATA_ROOT/dataset_metagenome_gut/fof.list:100:-fm"
    "tara:$DATA_ROOT/dataset_metagenome_tara/fof.list:10:-fq"
)

bench_init all_datasets "tuna kmc fastk"

for spec in "${DATASETS[@]}"; do
    IFS=: read -r ds fof maxn fmt <<< "$spec"
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    (( MAX_FILES > 0 && MAX_FILES < maxn )) && maxn=$MAX_FILES
    echo ""
    echo "== $ds  (m=${FORCE_M:-auto}, $fmt, first $maxn files) =="
    i=0
    used=()
    while IFS= read -r f && (( i < maxn )); do
        [[ -z "$f" || ! -f "$f" ]] && continue
        i=$((i+1))
        used+=("$f")
        echo "  [$i] $(basename "$f")"
        run_tuna "$ds" "$i" "$(basename "$f")" "$f"
        run_kmc  "$ds" "$i" "$(basename "$f")" "$f" "$fmt"
        # FastK takes a fof of extension-safe symlinks, so build a one-entry
        # one for this file rather than handing it the path directly.
        one="$AUX/onefof_${ds}_${i}.list"; printf '%s\n' "$f" > "$one"
        fkf="$AUX/fastkfof_${ds}_${i}.list"; fastk_fof "$one" 1 "$fkf"
        run_fastk "$ds" "$i" "$(basename "$f")" "$fkf"
        rm -f "$one" "$fkf"
    done < "$fof"

    # Same files again, this time as one collection.
    (( ${#used[@]} == 0 )) && continue
    sub="$AUX/subfof_${ds}.list"; printf '%s\n' "${used[@]}" > "$sub"
    echo ""
    echo "  [all] $ds as one collection (${#used[@]} files)"
    run_tuna "$ds" all "${#used[@]}" "@$sub"
    run_kmc  "$ds" all "${#used[@]}" "@$sub" "$fmt"
    fkf="$AUX/fastkfof_${ds}_all.list"; fastk_fof "$sub" "${#used[@]}" "$fkf"
    run_fastk "$ds" all "${#used[@]}" "$fkf"
    rm -f "$sub" "$fkf"
done

bench_done
