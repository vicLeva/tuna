#!/usr/bin/env bash
# bench_all_datasets.sh - per-file and whole-collection measurements.
#
# Two passes, which are two different experiments rather than two views of
# one:
#
#   per file     a sample of each collection counted one file at a time. The
#                "typical single input" view: what one genome, one metagenome
#                sample or one read set costs. A sample is enough, since the
#                figure reported is a median over the files.
#   whole        the collection counted in a single run, over as many files
#                as the paper claims for it. What a user pays to count the
#                whole thing at once, which is the regime where abundance is
#                high and peak memory matters.
#
# The two therefore have separate file lists and separate counts: a median
# over 100 genomes says the same thing as a median over 25000, but a
# collection run over 100 says nothing about one over 25000. Salmonella also
# draws from a different fof in each pass, the per-file sample from the small
# list and the collection from the pangenome one.
#
# Tools: tuna, KMC and FastK. FastK fails on a large share of these inputs
# (Tabex segfaults on the big ones), so expect fail rows rather than results
# on human and tara; the bin measurement still stands where only ascii fails.
#
# ONLY_PASS   "file", "all", or both when unset. Lets the two be run
#             separately, since the collection pass is far the more expensive
#             and is the one most likely to be repeated on its own.
#
# Output: $ROOT/all_datasets.csv
#   key = file index, label = file name   for the per-file rows
#   key = "all",      label = file count  for the whole-collection rows

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/bench_common.sh"

: "${ONLY_PASS:=}"
pass_enabled() { [[ -z "$ONLY_PASS" || " $ONLY_PASS " == *" $1 "* ]]; }

# name : fof : max files : kmc format flag
PERFILE=(
    "ecoli:$DATA_ROOT/dataset_genome_ecoli/fof.list:100:-fm"
    "human:$DATA_ROOT/dataset_genome_human/fof.list:10:-fm"
    "salmonella:$DATA_ROOT/dataset_pangenome_salmonella/fof.list:100:-fm"
    "gut:$DATA_ROOT/dataset_metagenome_gut/fof.list:100:-fm"
    "tara:$DATA_ROOT/dataset_metagenome_tara/fof.list:10:-fq"
)

# name : fof : files to take : kmc format flag
# The counts are what the paper reports for each collection. A fof holding
# fewer than that contributes everything it has.
COLLECTION=(
    "ecoli:$DATA_ROOT/dataset_genome_ecoli/fof.list:3500:-fm"
    "human:$DATA_ROOT/dataset_genome_human/fof.list:60:-fm"
    "salmonella:$DATA_ROOT/dataset_pangenome_salmonella/fof_150k.list:25000:-fm"
    "gut:$DATA_ROOT/dataset_metagenome_gut/fof.list:10000:-fm"
    "tara:$DATA_ROOT/dataset_metagenome_tara/fof.list:12:-fq"
)

bench_init all_datasets "tuna kmc fastk"

# --- pass 1: one file at a time --------------------------------------------
if pass_enabled file; then
for spec in "${PERFILE[@]}"; do
    IFS=: read -r ds fof maxn fmt <<< "$spec"
    dataset_enabled "$ds" || continue
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    (( MAX_FILES > 0 && MAX_FILES < maxn )) && maxn=$MAX_FILES
    echo ""
    echo "== $ds  (m=${FORCE_M:-auto}, $fmt, first $maxn files, one at a time) =="
    i=0
    while IFS= read -r f && (( i < maxn )); do
        [[ -z "$f" || ! -f "$f" ]] && continue
        i=$((i+1))
        echo "  [$i] $(basename "$f")"
        run_tuna "$ds" "$i" "$(basename "$f")" "$f"
        run_kmc  "$ds" "$i" "$(basename "$f")" "$f" "$fmt"
        # FastK takes a fof of extension-safe symlinks, so build a one-entry
        # one for this file rather than handing it the path directly. Skipped
        # when FastK is not running, so a tuna-only pass leaves no symlinks.
        if tool_enabled fastk; then
            one="$AUX/onefof_${ds}_${i}.list"; printf '%s\n' "$f" > "$one"
            fkf="$AUX/fastkfof_${ds}_${i}.list"; fastk_fof "$one" 1 "$fkf"
            run_fastk "$ds" "$i" "$(basename "$f")" "$fkf"
            rm -f "$one" "$fkf"
        fi
    done < "$fof"
done
fi

# --- pass 2: the whole collection in one run -------------------------------
if pass_enabled all; then
for spec in "${COLLECTION[@]}"; do
    IFS=: read -r ds fof want fmt <<< "$spec"
    dataset_enabled "$ds" || continue
    [[ -f "$fof" ]] || { echo "  [skip] $ds: no fof at $fof"; continue; }
    (( MAX_FILES > 0 && MAX_FILES < want )) && want=$MAX_FILES

    # Take the first $want paths that exist. A fof listing more than the
    # collection actually has on disk is common, so count what is there
    # rather than trusting the line count.
    sub="$AUX/subfof_${ds}.list"
    : > "$sub"
    n=0
    while IFS= read -r f && (( n < want )); do
        [[ -z "$f" || ! -f "$f" ]] && continue
        n=$((n+1)); printf '%s\n' "$f" >> "$sub"
    done < "$fof"
    (( n == 0 )) && { echo "  [skip] $ds: no readable file in $fof"; continue; }

    echo ""
    echo "== $ds  (m=${FORCE_M:-auto}, $fmt, $n files as one collection) =="
    (( n < want )) && echo "  [note] asked for $want, $fof holds $n"
    run_tuna "$ds" all "$n" "@$sub"
    run_kmc  "$ds" all "$n" "@$sub" "$fmt"
    if tool_enabled fastk; then
        fkf="$AUX/fastkfof_${ds}_all.list"; fastk_fof "$sub" "$n" "$fkf"
        run_fastk "$ds" all "$n" "$fkf"
        rm -f "$fkf"
    fi
    rm -f "$sub"
done
fi

bench_done
