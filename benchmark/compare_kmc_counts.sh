#!/usr/bin/env bash
#
# Benchmark Tuna and KMC without serializing k-mers, and require exact agreement
# on both total and distinct canonical k-mer counts.
#
# Usage:
#   compare_kmc_counts.sh INPUT [TUNA_OUTPUT_PLACEHOLDER]
#
# Configuration is supplied through environment variables:
#   TUNA=/path/to/tuna       KMC=/path/to/kmc
#   THREADS=8                RAM_GB=8
#   K=31                     M=15
#   CI=1                     CX=4294967295
#   KMC_FORMAT=-fq           WORK_ROOT=/scratch/path
#
# INPUT may be a FASTA/FASTQ path or a Tuna/KMC @file list. The second argument
# is only a CLI placeholder for Tuna in count-only mode and defaults to
# /dev/null.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 INPUT [TUNA_OUTPUT_PLACEHOLDER]" >&2
    exit 2
fi

INPUT=$1
TUNA_OUTPUT=${2:-/dev/null}
TUNA=${TUNA:-tuna}
KMC=${KMC:-kmc}
THREADS=${THREADS:-8}
RAM_GB=${RAM_GB:-8}
K=${K:-31}
M=${M:-15}
CI=${CI:-1}
CX=${CX:-4294967295}
KMC_FORMAT=${KMC_FORMAT:--fq}
WORK_ROOT=${WORK_ROOT:-${TMPDIR:-/tmp}}

TUNA=$(command -v "$TUNA" 2>/dev/null || printf '%s' "$TUNA")
KMC=$(command -v "$KMC" 2>/dev/null || printf '%s' "$KMC")
[[ -x "$TUNA" ]] || { echo "Tuna executable not found: $TUNA" >&2; exit 2; }
[[ -x "$KMC" ]] || { echo "KMC executable not found: $KMC" >&2; exit 2; }

run_dir=$(mktemp -d "$WORK_ROOT/tuna-kmc-compare.XXXXXX")
trap 'rm -rf "$run_dir"' EXIT

tuna_work="$run_dir/tuna"
kmc_work="$run_dir/kmc"
kmc_db="$run_dir/kmc-db"
tuna_log="$run_dir/tuna.log"
kmc_json="$run_dir/kmc.json"
mkdir -p "$tuna_work" "$kmc_work"

/usr/bin/time -f '%e' -o "$run_dir/tuna.wall" \
    "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" \
    -ci "$CI" -cx "$CX" -hp -co \
    -w "$tuna_work" "$INPUT" "$TUNA_OUTPUT" 2>"$tuna_log"

/usr/bin/time -f '%e' -o "$run_dir/kmc.wall" \
    "$KMC" -k"$K" -m"$RAM_GB" -t"$THREADS" -ci"$CI" -cx"$CX" -cs4294967295 \
    "$KMC_FORMAT" -hp -w -j"$kmc_json" \
    "$INPUT" "$kmc_db" "$kmc_work" >/dev/null 2>&1

value_after_colon() {
    local key=$1 file=$2
    awk -F: -v key="$key" '$1 == key {gsub(/[[:space:]s]/, "", $2); print $2; exit}' "$file"
}

json_integer() {
    local key=$1 file=$2
    awk -F: -v key="\"$key\"" '
        {
            lhs=$1
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", lhs)
            if (lhs == key) {
                value=$2
                gsub(/[^0-9]/, "", value)
                print value
                exit
            }
        }
    ' "$file"
}

tuna_total=$(value_after_colon total_kmers "$tuna_log")
tuna_unique=$(value_after_colon unique_kmers "$tuna_log")
tuna_p1=$(value_after_colon phase1 "$tuna_log")
tuna_p2=$(value_after_colon phase2 "$tuna_log")
tuna_wall=$(<"$run_dir/tuna.wall")

kmc_total=$(json_integer '#Total no. of k-mers' "$kmc_json")
kmc_unique=$(json_integer '#Unique_counted_k-mers' "$kmc_json")
kmc_p1=$(sed -n 's/.*"1st_stage":[[:space:]]*"\([^"]*\)".*/\1/p' "$kmc_json")
kmc_p2=$(sed -n 's/.*"2nd_stage":[[:space:]]*"\([^"]*\)".*/\1/p' "$kmc_json")
kmc_p1=${kmc_p1%s}
kmc_p2=${kmc_p2%s}
kmc_wall=$(<"$run_dir/kmc.wall")

for value_name in tuna_total tuna_unique kmc_total kmc_unique; do
    [[ -n ${!value_name} ]] || {
        echo "Failed to parse $value_name; retained logs are unavailable after exit." >&2
        exit 2
    }
done

printf '%-6s %10s %10s %20s %20s\n' tool wall_s phase1_s total_kmers distinct_kmers
printf '%-6s %10s %10s %20s %20s\n' tuna "$tuna_wall" "$tuna_p1" "$tuna_total" "$tuna_unique"
printf '%-6s %10s %10s %20s %20s\n' kmc "$kmc_wall" "$kmc_p1" "$kmc_total" "$kmc_unique"
printf 'phase2: tuna=%ss kmc=%ss\n' "$tuna_p2" "$kmc_p2"

if [[ "$tuna_total" != "$kmc_total" || "$tuna_unique" != "$kmc_unique" ]]; then
    echo "FAIL: Tuna and KMC counts differ." >&2
    exit 1
fi

echo "PASS: total and distinct counts match exactly."
