#!/usr/bin/env bash
# bench_k_sweep.sh — tuna vs KMC across k values on one human and one E.coli file
#
# Sweeps k = 31, 36, 41, ... up to K_MAX (step 5).
# For each k, runs m = 21 and m = k-10 (deduplicated when equal, e.g. k=31).
# One tuna binary is compiled per (k, m) pair with -DFIXED_K / -DFIXED_M and
# cached in BUILD_CACHE; use REBUILD=1 to force recompilation.
# KMC is run once per (dataset, k) since it does not depend on m.
#
# Usage: bash bench_k_sweep.sh [KEY=VALUE ...]
#
# Key=Value options:
#   THREADS=<int>   worker threads               [default: 4]
#   KMC_RAM=<int>   KMC RAM limit in GB          [default: 8]
#   K_MAX=<int>     highest k value to test      [default: 71]
#   REBUILD=1       force recompile all binaries [default: off]
#
# Path overrides (env vars):
#   TUNA_SRC=<dir>    tuna source directory     [default: directory above this script]
#   KMC=<path>        kmc binary                [default: kmc  (resolved via PATH)]
#   KMC_DUMP=<path>   kmc_dump binary           [default: kmc_dump (resolved via PATH)]
#   BUILD_CACHE=<dir> pre-built binary cache    [default: TUNA_SRC/build_k_sweep]
#   WORK=<dir>        scratch and results dir   [default: ~/tuna_bench_tmp]
#   HUMAN_FOF=<path>  human genome fof          [default: /WORKS/vlevallois/data/dataset_genome_human/fof.list]
#   ECOLI_FOF=<path>  E.coli genome fof         [default: /WORKS/vlevallois/data/dataset_genome_ecoli/fof.list]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parse KEY=VALUE arguments ─────────────────────────────────────────────────
for arg in "$@"; do
    case "$arg" in
        THREADS=*)  THREADS="${arg#THREADS=}" ;;
        KMC_RAM=*)  KMC_RAM="${arg#KMC_RAM=}" ;;
        K_MAX=*)    K_MAX="${arg#K_MAX=}" ;;
        REBUILD=*)  REBUILD="${arg#REBUILD=}" ;;
        *) echo "bench_k_sweep: unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS="${THREADS:-4}"
KMC_RAM="${KMC_RAM:-8}"
K_MAX="${K_MAX:-71}"
REBUILD="${REBUILD:-0}"

TUNA_SRC="${TUNA_SRC:-$SCRIPT_DIR/..}"
KMC="${KMC:-kmc}"
KMC_DUMP="${KMC_DUMP:-kmc_dump}"
BUILD_CACHE="${BUILD_CACHE:-$TUNA_SRC/build_k_sweep}"
WORK="${WORK:-$HOME/tuna_bench_tmp}"
HUMAN_FOF="${HUMAN_FOF:-/WORKS/vlevallois/data/dataset_genome_human/fof.list}"
ECOLI_FOF="${ECOLI_FOF:-/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list}"

RESULTS="$WORK/bench_k_sweep_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS" "$BUILD_CACHE"

# ── Validate external tools ───────────────────────────────────────────────────
KMC=$(command -v "$KMC"      2>/dev/null || echo "$KMC")
KMC_DUMP=$(command -v "$KMC_DUMP" 2>/dev/null || echo "$KMC_DUMP")
[[ -x "$KMC"      ]] || { echo "ERROR: kmc not found: $KMC"      >&2; exit 1; }
[[ -x "$KMC_DUMP" ]] || { echo "ERROR: kmc_dump not found: $KMC_DUMP" >&2; exit 1; }
/usr/bin/time -v true 2>/dev/null \
    || { echo "ERROR: /usr/bin/time -v not available (need GNU time)" >&2; exit 1; }

# ── Build (k, m) sweep table ──────────────────────────────────────────────────
# K values: 31, 36, 41, ..., K_MAX (step 5).
# M values per k: 21 and k-10, deduplicated.
K_VALUES=()
KM_PAIRS=()   # "k:m" strings for all unique pairs
declare -A K_M_LIST  # k -> space-separated list of m values

k=31
while (( k <= K_MAX )); do
    K_VALUES+=("$k")
    ms=""
    for candidate in 17 21 25; do
        (( candidate < k )) && ms="$ms $candidate"
    done
    ms="${ms# }"  # trim leading space
    K_M_LIST[$k]="$ms"
    for m in $ms; do KM_PAIRS+=("$k:$m"); done
    (( k += 5 ))
done

# ── Dataset registry ──────────────────────────────────────────────────────────
# Format: "name:fof_path:kmc_format"
DATASETS=(
    "ecoli:$ECOLI_FOF:-fm"
    "human:$HUMAN_FOF:-fm"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

wall_to_s() {
    awk -F: '{if(NF==3) printf "%.2f", $1*3600+$2*60+$3;
              else       printf "%.2f", $1*60+$2}' <<< "$1"
}

wall_from_file() {
    local t
    t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    wall_to_s "$t"
}

rss_mb() {
    local kb
    kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

HAS_SUDO=false
sudo -n true 2>/dev/null && HAS_SUDO=true

drop_caches() {
    if $HAS_SUDO; then
        sync
        sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
    fi
}

# ── Phase 0: build tuna binaries ──────────────────────────────────────────────
echo "=== Building tuna binaries (FIXED_K + FIXED_M) ==="
echo "    source: $TUNA_SRC"
echo "    cache:  $BUILD_CACHE"
echo ""

declare -A TUNA_BIN  # "k:m" -> binary path

for pair in "${KM_PAIRS[@]}"; do
    k="${pair%%:*}"
    m="${pair##*:}"
    build_dir="$BUILD_CACHE/k${k}_m${m}"
    binary="$build_dir/tuna"

    if [[ -x "$binary" && "$REBUILD" != "1" ]]; then
        echo "  [cached]  k=$k m=$m  ($binary)"
        TUNA_BIN["$pair"]="$binary"
        continue
    fi

    echo "  [build]   k=$k m=$m ..."
    mkdir -p "$build_dir"
    cmake -S "$TUNA_SRC" -B "$build_dir" \
        -DFIXED_K="$k" -DFIXED_M="$m" \
        -DCMAKE_BUILD_TYPE=Release \
        > "$build_dir/cmake.log" 2>&1 \
        || { echo "    cmake FAILED — check $build_dir/cmake.log" >&2; exit 1; }
    cmake --build "$build_dir" --target tuna -j"$(nproc)" \
        > "$build_dir/make.log" 2>&1 \
        || { echo "    make FAILED — check $build_dir/make.log" >&2; exit 1; }
    echo "    done → $binary"
    TUNA_BIN["$pair"]="$binary"
done

echo ""

# ── CSV header ────────────────────────────────────────────────────────────────
CSV="$RESULTS/bench_k_sweep.csv"
echo "dataset,filename,k,m,tool,wall_s,rss_mb,unique_kmers,phase1_s,phase2_s" > "$CSV"

# ── Phase 1: run benchmarks ───────────────────────────────────────────────────
echo "=== Running benchmarks ==="
echo "    threads=$THREADS  kmc_ram=${KMC_RAM}GB"
echo "    k values: ${K_VALUES[*]}"
if $HAS_SUDO; then
    echo "    cache drops: ENABLED"
else
    echo "    cache drops: DISABLED (no passwordless sudo)"
fi
echo ""

for ENTRY in "${DATASETS[@]}"; do
    IFS=: read -r DS_NAME FOF KMC_FMT <<< "$ENTRY"

    if [[ ! -f "$FOF" ]]; then
        echo "  [SKIP] $DS_NAME: fof not found: $FOF"
        continue
    fi

    FILE=$(head -1 "$FOF")
    if [[ -z "$FILE" || ! -f "$FILE" ]]; then
        echo "  [SKIP] $DS_NAME: first file from fof not found: $FILE"
        continue
    fi
    FNAME=$(basename "$FILE")

    echo "──── Dataset: $DS_NAME  file: $FNAME ────"

    DS_DIR="$RESULTS/$DS_NAME"
    mkdir -p "$DS_DIR"

    for k in "${K_VALUES[@]}"; do
        echo "  k=$k"

        # ── KMC (run once per k, m-independent) ───────────────────────────────
        kmc_db="$WORK/kmc_db_${DS_NAME}_k${k}_$$"
        kmc_tmp="$WORK/kmc_tmp_${DS_NAME}_k${k}_$$"
        kmc_dump_out=$(mktemp "$WORK/kmc_dump_XXXXXX.tsv")
        kmc_timefile="$DS_DIR/kmc_k${k}.timefile"
        kmc_dump_timefile="$DS_DIR/kmc_k${k}_dump.timefile"
        mkdir -p "$kmc_tmp"

        drop_caches
        /usr/bin/time -v -o "$kmc_timefile" \
            "$KMC" -k"$k" -m"$KMC_RAM" -ci1 "$KMC_FMT" -t"$THREADS" \
            "$FILE" "$kmc_db" "$kmc_tmp" \
            > "$DS_DIR/kmc_k${k}.log" 2>&1 \
            || { echo "    [kmc k=$k] FAILED — check $DS_DIR/kmc_k${k}.log"; rm -rf "$kmc_tmp"; continue; }
        rm -rf "$kmc_tmp"

        drop_caches
        /usr/bin/time -v -o "$kmc_dump_timefile" \
            "$KMC_DUMP" -ci1 "$kmc_db" "$kmc_dump_out" \
            >> "$DS_DIR/kmc_k${k}.log" 2>&1 \
            || { echo "    [kmc_dump k=$k] FAILED"; rm -f "${kmc_db}.kmc_pre" "${kmc_db}.kmc_suf" "$kmc_dump_out"; continue; }
        rm -f "${kmc_db}.kmc_pre" "${kmc_db}.kmc_suf"

        k_count_wall=$(wall_from_file "$kmc_timefile")
        k_dump_wall=$(wall_from_file  "$kmc_dump_timefile")
        k_wall=$(awk "BEGIN{printf \"%.2f\", $k_count_wall + $k_dump_wall}")
        k_rss=$(rss_mb "$kmc_timefile")
        K_KMERS=$(wc -l < "$kmc_dump_out" || echo 0)
        rm -f "$kmc_dump_out"

        echo "    [kmc]       total=${k_wall}s (count=${k_count_wall}s dump=${k_dump_wall}s)  RSS=${k_rss}MB  kmers=${K_KMERS}"
        echo "${DS_NAME},${FNAME},${k},na,kmc,${k_wall},${k_rss},${K_KMERS},${k_count_wall},${k_dump_wall}" >> "$CSV"

        # ── tuna (one run per m value for this k) ─────────────────────────────
        IFS=' ' read -ra M_VALUES <<< "${K_M_LIST[$k]}"
        for m in "${M_VALUES[@]}"; do
            pair="${k}:${m}"
            tuna_bin="${TUNA_BIN[$pair]}"
            tuna_out=$(mktemp "$WORK/tuna_out_XXXXXX.tsv")
            tuna_work="$WORK/tuna_work_${DS_NAME}_k${k}_m${m}_$$"
            tuna_timefile="$DS_DIR/tuna_k${k}_m${m}.timefile"
            tuna_stderr="$DS_DIR/tuna_k${k}_m${m}.stderr"
            mkdir -p "$tuna_work"

            drop_caches
            /usr/bin/time -v -o "$tuna_timefile" \
                "$tuna_bin" -k "$k" -m "$m" -t "$THREADS" \
                -w "$tuna_work/" "$FILE" "$tuna_out" \
                2>"$tuna_stderr" \
                || { echo "    [tuna k=$k m=$m] FAILED — check $tuna_stderr"; rm -rf "$tuna_work" "$tuna_out"; continue; }
            rm -rf "$tuna_work"

            t_wall=$(wall_from_file "$tuna_timefile")
            t_rss=$(rss_mb "$tuna_timefile")
            t_kmers=$(wc -l < "$tuna_out" || echo 0)
            t_p1=$(grep "^phase1:" "$tuna_stderr" | awk -F: '{gsub(/s/,"",$2); printf "%.3f",$2}' || echo na)
            t_p2=$(grep "^phase2:" "$tuna_stderr" | awk -F: '{gsub(/s/,"",$2); printf "%.3f",$2}' || echo na)
            rm -f "$tuna_out"

            # Correctness check vs KMC
            if [[ "$t_kmers" -eq "$K_KMERS" ]] 2>/dev/null; then
                match="OK"
            else
                match="DIFF(tuna=$t_kmers kmc=$K_KMERS)"
            fi

            echo "    [tuna m=$m]  total=${t_wall}s (p1=${t_p1}s p2=${t_p2}s)  RSS=${t_rss}MB  kmers=${t_kmers}  [$match]"
            echo "${DS_NAME},${FNAME},${k},${m},tuna,${t_wall},${t_rss},${t_kmers},${t_p1},${t_p2}" >> "$CSV"
        done

        echo ""
    done
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo "=== Done ==="
echo "Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
echo ""
column -t -s, "$CSV"
