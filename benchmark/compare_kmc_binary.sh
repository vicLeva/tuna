#!/usr/bin/env bash
#
# Benchmark real binary count output and require exact agreement on total and
# distinct canonical k-mer counts. Tuna writes KFF; KMC writes its native DB.
#
# Usage:
#   compare_kmc_binary.sh INPUT RESULT_DIR
#
# Environment:
#   TUNA=/path/to/tuna       KMC=/path/to/kmc
#   THREADS=16               RAM_GB=256
#   K=31                     M=15
#   CI=1                     CX=4294967295
#   KMC_FORMAT=-fq

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 INPUT RESULT_DIR" >&2
    exit 2
fi

INPUT=$1
RESULT_DIR=$2
TUNA=${TUNA:-tuna}
KMC=${KMC:-kmc}
THREADS=${THREADS:-16}
RAM_GB=${RAM_GB:-256}
K=${K:-31}
M=${M:-15}
CI=${CI:-1}
CX=${CX:-4294967295}
KMC_FORMAT=${KMC_FORMAT:--fq}

TUNA=$(command -v "$TUNA" 2>/dev/null || printf '%s' "$TUNA")
KMC=$(command -v "$KMC" 2>/dev/null || printf '%s' "$KMC")
[[ -x "$TUNA" ]] || { echo "Tuna executable not found: $TUNA" >&2; exit 2; }
[[ -x "$KMC" ]] || { echo "KMC executable not found: $KMC" >&2; exit 2; }

if [[ -e "$RESULT_DIR" ]] && [[ -n $(find "$RESULT_DIR" -mindepth 1 -print -quit) ]]; then
    echo "Result directory is not empty: $RESULT_DIR" >&2
    exit 2
fi

mkdir -p "$RESULT_DIR/tuna_work" "$RESULT_DIR/kmc_work"
tuna_out="$RESULT_DIR/tuna.kff"
kmc_db="$RESULT_DIR/kmc"

/usr/bin/time -f '%e' -o "$RESULT_DIR/tuna.wall" \
    "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" \
    -ci "$CI" -cx "$CX" -hp -kff -w "$RESULT_DIR/tuna_work" \
    "$INPUT" "$tuna_out" 2>"$RESULT_DIR/tuna.log"

/usr/bin/time -f '%e' -o "$RESULT_DIR/kmc.wall" \
    "$KMC" -k"$K" -m"$RAM_GB" -t"$THREADS" -ci"$CI" -cx"$CX" \
    -cs4294967295 "$KMC_FORMAT" -hp -j"$RESULT_DIR/kmc.json" \
    "$INPUT" "$kmc_db" "$RESULT_DIR/kmc_work" \
    >"$RESULT_DIR/kmc.stdout" 2>"$RESULT_DIR/kmc.stderr"

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

tuna_total=$(value_after_colon total_kmers "$RESULT_DIR/tuna.log")
tuna_unique=$(value_after_colon unique_kmers "$RESULT_DIR/tuna.log")
kmc_total=$(json_integer '#Total no. of k-mers' "$RESULT_DIR/kmc.json")
kmc_unique=$(json_integer '#Unique_counted_k-mers' "$RESULT_DIR/kmc.json")

if [[ "$tuna_total" != "$kmc_total" || "$tuna_unique" != "$kmc_unique" ]]; then
    echo "FAIL: Tuna and KMC counts differ." >&2
    exit 1
fi

tuna_bytes=$(stat -c %s "$tuna_out")
kmc_bytes=$((
    $(stat -c %s "${kmc_db}.kmc_pre") +
    $(stat -c %s "${kmc_db}.kmc_suf")
))

printf '%-6s %10s %20s %20s %16s\n' tool wall_s total_kmers distinct_kmers output_bytes
printf '%-6s %10s %20s %20s %16s\n' \
    tuna "$(<"$RESULT_DIR/tuna.wall")" "$tuna_total" "$tuna_unique" "$tuna_bytes"
printf '%-6s %10s %20s %20s %16s\n' \
    kmc "$(<"$RESULT_DIR/kmc.wall")" "$kmc_total" "$kmc_unique" "$kmc_bytes"
echo "PASS: exact summary counts match with real binary output."
