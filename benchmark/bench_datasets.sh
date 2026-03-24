#!/usr/bin/env bash
# bench_datasets.sh — tuna (m=21) vs KMC on 6 datasets
#
# Usage: bash bench_datasets.sh [THREADS] [K] [KMC_RAM_GB] [KMC_CACHE_CSV]
# Example: bash bench_datasets.sh 8 31 250 /path/to/old/bench.csv
#
# KMC_CACHE_CSV: path to a previous bench.csv that already contains KMC
#   per-file rows.  When set, run_kmc will look up (dataset, filename) in
#   that file and copy the cached row instead of re-running KMC.  KMC is
#   always executed for the full-dataset experiment (fname="full").
#
# Two experiments per dataset:
#
#   Per-file   — N files selected from the fof (head or spread mode), one
#                tuna/KMC call per file.  Measures single-file throughput.
#
#   Full       — One tuna/KMC call on the complete fof.list (all files at
#                once).  Measures multi-file parallel throughput.
#                tuna receives "@fof.list"; KMC likewise.
#
# Sampling modes (per-file only):
#   head   — take the first N files in the fof (default)
#   spread — pick N files evenly spread across the full fof (always
#             includes first and last)
#
# Output layout:
#   $RESULTS/bench.csv            — one row per (dataset, file, tool, m)
#   $RESULTS/<ds>/<tag>.stderr    — tuna structured timing lines
#   $RESULTS/<ds>/<tag>.timefile  — GNU time -v output
#   $RESULTS/<ds>/<tag>.kmc.log   — KMC stdout/stderr
#
# file_idx=0 / filename="full" marks the whole-dataset rows in the CSV.
#
# k-mer output files are written to a temp path and deleted after counting
# to avoid accumulating hundreds of GB of TSV.

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────

THREADS=${1:-8}
K=${2:-31}
KMC_RAM=${3:-250}      # GB — KMC -m flag
KMC_CACHE_CSV=${4:-}   # optional: path to prior bench.csv with KMC per-file rows

M_VALUES=(21)

TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
KMC=kmc
KMC_DUMP=kmc_dump
WORK=/WORKS/vlevallois/test_tuna
RESULTS="$WORK/bench_datasets_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$RESULTS"

# ── Dataset registry ──────────────────────────────────────────────────────────
# Format: "name:fof_path:kmc_format:max_files:mode:full_max"
#   max_files — how many files to use for the per-file experiment
#   mode      — head (default) or spread
#   full_max  — how many files for the full-dataset experiment (0 = all)
DATASETS=(
    "ecoli:/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list:-fm:100:head:3000"
    "human:/WORKS/vlevallois/data/dataset_genome_human/fof.list:-fm:10:head:30"
    "salmonella:/WORKS/vlevallois/data/dataset_pangenome_salmonella/fof.list:-fm:100:head:3000"
    "gut:/WORKS/vlevallois/data/dataset_metagenome_gut/fof.list:-fm:100:head:3000"
    "tara:/WORKS/vlevallois/data/dataset_metagenome_tara/fof.list:-fq:10:head:0"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

# GNU time "Elapsed" field (M:SS.ss or H:MM:SS.ss) → seconds
wall_to_s() {
    awk -F: '{if(NF==3) printf "%.2f", $1*3600+$2*60+$3;
              else       printf "%.2f", $1*60+$2}' <<< "$1"
}

rss_mb() {
    local kb
    kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

wall_from_file() {
    local t
    t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    wall_to_s "$t"
}

# ── CSV header ────────────────────────────────────────────────────────────────

CSV="$RESULTS/bench.csv"
echo "dataset,file_idx,filename,tool,m,wall_s,rss_mb,unique_kmers,phase1_s,phase2_s" > "$CSV"
# For tuna: phase1_s=partitioning, phase2_s=counting+writing
# For kmc:  phase1_s=kmc_count_wall, phase2_s=kmc_dump_wall, wall_s=sum
# file_idx=0 / filename=full → whole-dataset row

echo "=== Dataset benchmark ==="
echo "    k=$K  m=${M_VALUES[*]}  threads=$THREADS  kmc_ram=${KMC_RAM}GB"
echo "    results: $RESULTS"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────

for ENTRY in "${DATASETS[@]}"; do

    IFS=: read -r DS_NAME FOF KMC_FMT DS_MAX DS_MODE DS_FULL_MAX <<< "$ENTRY"
    DS_MAX="${DS_MAX:-10}"
    DS_MODE="${DS_MODE:-head}"
    DS_FULL_MAX="${DS_FULL_MAX:-0}"

    if [ ! -f "$FOF" ]; then
        echo "  [SKIP] $DS_NAME: fof not found: $FOF"
        continue
    fi

    TOTAL=$(wc -l < "$FOF")
    N=$(( TOTAL < DS_MAX ? TOTAL : DS_MAX ))

    DS_DIR="$RESULTS/$DS_NAME"
    mkdir -p "$DS_DIR"

    echo "──── Dataset: $DS_NAME  ($N / $TOTAL files, format: $KMC_FMT, mode: $DS_MODE) ────────"

    if [ "$DS_MODE" = "spread" ] && [ "$N" -ge 2 ] && [ "$TOTAL" -gt "$N" ]; then
        mapfile -t ALL_FILES < "$FOF"
        FILES=()
        for i in $(seq 0 $(( N - 1 ))); do
            idx=$(( i * (TOTAL - 1) / (N - 1) ))
            FILES+=("${ALL_FILES[$idx]}")
        done
        unset ALL_FILES
    else
        mapfile -t FILES < <(head -n "$N" "$FOF")
    fi

    # Defined here so they share DS_DIR, KMC_FMT, CSV, etc. via closure.

    run_kmc() {
        local tag="$1" file="$2" fidx="$3" fname="$4"

        # ── Cache lookup (per-file only) ───────────────────────────────────────
        # For full-dataset rows (fname="full") always run KMC — no prior data.
        if [ -n "$KMC_CACHE_CSV" ] && [ -f "$KMC_CACHE_CSV" ] && [ "$fname" != "full" ]; then
            local cached
            cached=$(awk -F, -v ds="$DS_NAME" -v fn="$fname" \
                '$1==ds && $3==fn && $4=="kmc" {print; exit}' "$KMC_CACHE_CSV")
            if [ -n "$cached" ]; then
                # Columns: dataset,file_idx,filename,tool,m,wall_s,rss_mb,unique_kmers,phase1_s,phase2_s
                local k_wall k_rss k_p1 k_p2
                k_wall=$(echo "$cached" | cut -d, -f6)
                k_rss=$(echo  "$cached" | cut -d, -f7)
                K_KMERS=$(echo "$cached" | cut -d, -f8)
                k_p1=$(echo   "$cached" | cut -d, -f9)
                k_p2=$(echo   "$cached" | cut -d, -f10)
                echo "    [kmc]  total=${k_wall}s (count=${k_p1}s dump=${k_p2}s)  RSS=${k_rss}MB  kmers=${K_KMERS}  [cached]"
                echo "${DS_NAME},${fidx},${fname},kmc,na,${k_wall},${k_rss},${K_KMERS},${k_p1},${k_p2}" >> "$CSV"
                return
            fi
            echo "    [kmc]  no cache hit for $DS_NAME/$fname — running KMC"
        fi

        # ── Run KMC ───────────────────────────────────────────────────────────
        local kmc_db="$WORK/kmc_db_${tag}"
        local kmc_tmp="$WORK/kmc_tmp_${tag}"
        local kmc_dump_out
        kmc_dump_out=$(mktemp "$WORK/kmc_dump_XXXXXX.tsv")
        mkdir -p "$kmc_tmp"

        /usr/bin/time -v -o "$DS_DIR/${tag}.kmc.timefile" \
            "$KMC" -k"$K" -m"$KMC_RAM" -ci1 "$KMC_FMT" -t"$THREADS" \
            "$file" "$kmc_db" "$kmc_tmp" \
            > "$DS_DIR/${tag}.kmc.log" 2>&1 || {
                echo "    [kmc] FAILED — check $DS_DIR/${tag}.kmc.log"
            }
        rm -rf "$kmc_tmp"

        /usr/bin/time -v -o "$DS_DIR/${tag}.kmc_dump.timefile" \
            "$KMC_DUMP" -ci1 "$kmc_db" "$kmc_dump_out" \
            >> "$DS_DIR/${tag}.kmc.log" 2>&1 || {
                echo "    [kmc_dump] FAILED — check $DS_DIR/${tag}.kmc.log"
            }
        rm -f "${kmc_db}.kmc_pre" "${kmc_db}.kmc_suf"

        local k_count_wall k_dump_wall k_wall k_rss
        k_count_wall=$(wall_from_file "$DS_DIR/${tag}.kmc.timefile")
        k_dump_wall=$(wall_from_file "$DS_DIR/${tag}.kmc_dump.timefile")
        k_wall=$(awk "BEGIN{printf \"%.2f\", $k_count_wall + $k_dump_wall}")
        k_rss=$(rss_mb "$DS_DIR/${tag}.kmc.timefile")
        K_KMERS=$([ -f "$kmc_dump_out" ] && wc -l < "$kmc_dump_out" || echo 0)
        rm -f "$kmc_dump_out"

        echo "    [kmc]  total=${k_wall}s (count=${k_count_wall}s dump=${k_dump_wall}s)  RSS=${k_rss}MB  kmers=${K_KMERS}"
        echo "${DS_NAME},${fidx},${fname},kmc,na,${k_wall},${k_rss},${K_KMERS},${k_count_wall},${k_dump_wall}" >> "$CSV"
    }

    run_tuna() {
        local tag="$1" file="$2" fidx="$3" fname="$4"
        for M in "${M_VALUES[@]}"; do
            local tuna_out tuna_work
            tuna_out=$(mktemp "$WORK/tuna_out_XXXXXX.tsv")
            tuna_work="$WORK/tuna_work_${tag}_m${M}"
            mkdir -p "$tuna_work"

            /usr/bin/time -v -o "$DS_DIR/${tag}.m${M}.timefile" \
                "$TUNA" -k "$K" -m "$M" -t "$THREADS" \
                -w "$tuna_work/" "$file" "$tuna_out" \
                2>"$DS_DIR/${tag}.m${M}.stderr" || {
                    echo "    [tuna m=$M] FAILED — check $DS_DIR/${tag}.m${M}.stderr"
                }
            rm -rf "$tuna_work"

            local t_wall t_rss t_kmers t_p1 t_p2
            t_wall=$(wall_from_file "$DS_DIR/${tag}.m${M}.timefile")
            t_rss=$(rss_mb "$DS_DIR/${tag}.m${M}.timefile")
            t_kmers=$([ -f "$tuna_out" ] && wc -l < "$tuna_out" || echo 0)
            t_p1=$(grep "^phase1:" "$DS_DIR/${tag}.m${M}.stderr" \
                   | awk -F: '{gsub(/s/,"",$2); printf "%.3f", $2}' || echo na)
            t_p2=$(grep "^phase2:" "$DS_DIR/${tag}.m${M}.stderr" \
                   | awk -F: '{gsub(/s/,"",$2); printf "%.3f", $2}' || echo na)
            rm -f "$tuna_out"

            echo "    [tuna m=$M]  wall=${t_wall}s (p1=${t_p1}s p2=${t_p2}s)  RSS=${t_rss}MB  kmers=${t_kmers}"
            echo "${DS_NAME},${fidx},${fname},tuna,${M},${t_wall},${t_rss},${t_kmers},${t_p1},${t_p2}" >> "$CSV"

            TUNA_KMERS[$M]=$t_kmers
        done
    }

    # ── Per-file experiment ────────────────────────────────────────────────────

    for IDX in "${!FILES[@]}"; do

        FILE="${FILES[$IDX]}"
        FNAME=$(basename "$FILE")
        FIDX=$(( IDX + 1 ))
        TAG="${DS_NAME}_f$(printf '%04d' $FIDX)"
        K_KMERS=0        # set by run_kmc
        declare -A TUNA_KMERS  # set by run_tuna: TUNA_KMERS[m]=count

        # Alternate order each file so neither tool consistently gets warm cache.
        if (( FIDX % 2 == 1 )); then
            echo "  [$FIDX/$N] $FNAME  (kmc→tuna)"
            run_kmc  "$TAG" "$FILE" "$FIDX" "$FNAME"
            run_tuna "$TAG" "$FILE" "$FIDX" "$FNAME"
        else
            echo "  [$FIDX/$N] $FNAME  (tuna→kmc)"
            run_tuna "$TAG" "$FILE" "$FIDX" "$FNAME"
            run_kmc  "$TAG" "$FILE" "$FIDX" "$FNAME"
        fi

        # Correctness check — only print on mismatch.
        any_diff=0
        for M in "${M_VALUES[@]}"; do
            t_kmers=${TUNA_KMERS[$M]:-0}
            if [ "$t_kmers" -ne "$K_KMERS" ] 2>/dev/null; then
                echo "    [DIFF] tuna m=$M: tuna=$t_kmers  kmc=$K_KMERS"
                any_diff=1
            fi
        done
        [ $any_diff -eq 0 ] && echo "    [OK]  all m values match kmc ($K_KMERS unique k-mers)"

        unset TUNA_KMERS
    done

    # ── Full-dataset experiment ────────────────────────────────────────────────
    # One tuna call + one KMC call on (up to DS_FULL_MAX, or all) files.
    # tuna and KMC both accept "@fof.list" as multi-file input.

    if [ "$DS_FULL_MAX" -gt 0 ] && [ "$TOTAL" -gt "$DS_FULL_MAX" ]; then
        FULL_N=$DS_FULL_MAX
        FULL_FOF=$(mktemp "$WORK/full_fof_${DS_NAME}_XXXXXX.list")
        head -n "$DS_FULL_MAX" "$FOF" > "$FULL_FOF"
    else
        FULL_N=$TOTAL
        FULL_FOF="$FOF"
    fi

    echo ""
    echo "  [full] $DS_NAME — $FULL_N files"
    FTAG="${DS_NAME}_full"
    K_KMERS=0
    declare -A TUNA_KMERS

    run_kmc  "$FTAG" "@${FULL_FOF}" 0 "full"
    run_tuna "$FTAG" "@${FULL_FOF}" 0 "full"

    [ "$FULL_FOF" != "$FOF" ] && rm -f "$FULL_FOF"

    any_diff=0
    for M in "${M_VALUES[@]}"; do
        t_kmers=${TUNA_KMERS[$M]:-0}
        if [ "$t_kmers" -ne "$K_KMERS" ] 2>/dev/null; then
            echo "    [DIFF] tuna m=$M: tuna=$t_kmers  kmc=$K_KMERS"
            any_diff=1
        fi
    done
    [ $any_diff -eq 0 ] && echo "    [OK]  all m values match kmc ($K_KMERS unique k-mers)"

    unset TUNA_KMERS
    echo ""
done

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=== Done ==="
echo "Results: $RESULTS/bench.csv  ($(wc -l < "$CSV") rows)"
echo ""
column -t -s, "$CSV" | head -40
echo "..."
