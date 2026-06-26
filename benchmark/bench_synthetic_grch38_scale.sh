#!/usr/bin/env bash
# Synthetic GRCh38-derived Tuna vs KMC scaling benchmark.
#
# Generates one medium FASTA from a real human reference, then benchmarks:
#   repeated-files: N distinct symlinked input files with identical content
#   concat:         one FASTA containing N copies of the same content
#
# Timing reported:
#   core: Tuna to /dev/null; KMC count only
#   e2e:  Tuna TSV to scratch; KMC count + kmc_dump TSV to scratch
#
# All generated data, outputs, logs, and temporary files stay under WORK.

set -euo pipefail

REF=${REF:-/scratch2/tmp/grch38_1019/genome.fa.gz}
WORK=${WORK:-/scratch3/tmp/tuna-optimization-codex/bench/synthetic_grch38}
TMPDIR=${TMPDIR:-/scratch3/tmp/tuna-optimization-codex/tmp}

TUNA=${TUNA:-/scratch3/tmp/tuna-optimization-codex/build/tuna-k31-m21/tuna}
KMC=${KMC:-/scratch3/rob/tuna-optimization/KMC/bin/kmc}
KMC_DUMP=${KMC_DUMP:-/scratch3/rob/tuna-optimization/KMC/bin/kmc_dump}

K=${K:-31}
M=${M:-21}
THREADS=${THREADS:-8}
PARTS=${PARTS:-}
KMC_RAM_GB=${KMC_RAM_GB:-64}
BASE_MB=${BASE_MB:-64}
N_VALUES=${N_VALUES:-"1 2 4 8 16 32 64"}
SCENARIOS=${SCENARIOS:-"repeated concat"}
REPS=${REPS:-1}

mkdir -p "$WORK" "$TMPDIR"
export TMPDIR

DATA_DIR="$WORK/data"
RUN_DIR="$WORK/runs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$DATA_DIR" "$RUN_DIR"

BASE_FA="$DATA_DIR/grch38_${BASE_MB}mb.fa"
CSV="$RUN_DIR/bench_synthetic_grch38.csv"

require_exe() {
    local p="$1"
    [[ -x "$p" ]] || { echo "ERROR: not executable: $p" >&2; exit 1; }
}

require_exe "$TUNA"
require_exe "$KMC"
require_exe "$KMC_DUMP"
[[ -f "$REF" ]] || { echo "ERROR: missing reference: $REF" >&2; exit 1; }

echo "tool,scenario,n,rep,core_wall_s,core_rss_mb,e2e_wall_s,e2e_rss_mb,phase1_s,phase2_s,unique_kmers,n_parts,superkmers,input_bytes,notes" > "$CSV"

wall_to_s() {
    awk -F: '{if (NF == 3) printf "%.3f", $1*3600+$2*60+$3; else printf "%.3f", $1*60+$2}'
}

wall_from_time() {
    grep "Elapsed (wall clock)" "$1" | awk '{print $NF}' | wall_to_s
}

rss_mb_from_time() {
    local kb
    kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}

parse_tuna_key() {
    local f="$1" key="$2"
    grep "^${key}:" "$f" | awk -F': ' '{gsub(/s$/,"",$2); print $2}' | head -1
}

sum_s() {
    awk "BEGIN{printf \"%.3f\", $1 + $2}"
}

max_i() {
    awk "BEGIN{print ($1 > $2) ? $1 : $2}"
}

generate_base() {
    if [[ -s "$BASE_FA" ]]; then
        return
    fi
    echo "[setup] generating $BASE_FA from $REF (${BASE_MB} MiB ACTG)"
    python3 - "$REF" "$BASE_FA" "$BASE_MB" <<'PY'
import gzip
import sys

ref, out, mb = sys.argv[1], sys.argv[2], int(sys.argv[3])
target = mb * 1024 * 1024
written = 0
width = 80
buf = []

def flush(f):
    global buf, written
    if not buf:
        return
    s = ''.join(buf)
    for i in range(0, len(s), width):
        f.write(s[i:i+width] + '\n')
    written += len(s)
    buf = []

with gzip.open(ref, "rt", encoding="ascii", errors="ignore") as inp, open(out, "w") as out_f:
    out_f.write(f">grch38_first_actg_{mb}mb\n")
    for line in inp:
        if written >= target:
            break
        if line.startswith(">"):
            continue
        seq = ''.join(ch for ch in line.upper() if ch in "ACGT")
        if not seq:
            continue
        need = target - written - sum(map(len, buf))
        if need <= 0:
            flush(out_f)
            break
        buf.append(seq[:need])
        if sum(map(len, buf)) >= 1 << 20:
            flush(out_f)
    flush(out_f)

if written < target:
    raise SystemExit(f"reference ended after {written} bases, wanted {target}")
PY
}

prepare_repeated_fof() {
    local n="$1"
    local dir="$DATA_DIR/repeated_n${n}"
    local fof="$DATA_DIR/repeated_n${n}.fof"
    mkdir -p "$dir"
    : > "$fof"
    for i in $(seq 1 "$n"); do
        local p="$dir/copy_$(printf '%05d' "$i").fa"
        ln -sf "$BASE_FA" "$p"
        echo "$p" >> "$fof"
    done
    echo "$fof"
}

prepare_concat_file() {
    local n="$1"
    local out="$DATA_DIR/concat_n${n}.fa"
    if [[ ! -s "$out" ]]; then
        echo "[setup] creating concat input $out" >&2
        : > "$out"
        for i in $(seq 1 "$n"); do
            cat "$BASE_FA" >> "$out"
        done
    fi
    echo "$out"
}

run_tuna() {
    local tag="$1" input_arg="$2" out_path="$3" time_f="$4" stderr_f="$5"
    local work_d="$RUN_DIR/${tag}.tuna_work"
    mkdir -p "$work_d"
    local part_args=()
    local count_args=()
    if [[ -n "$PARTS" ]]; then
        part_args=(-n "$PARTS")
    fi
    if [[ "$out_path" == "/dev/null" ]]; then
        count_args=(-co)
    fi
    /usr/bin/time -v -o "$time_f" \
        "$TUNA" -k "$K" -m "$M" -t "$THREADS" "${part_args[@]}" "${count_args[@]}" -hp -w "$work_d" \
        "$input_arg" "$out_path" \
        2>"$stderr_f"
    rm -rf "$work_d"
}

run_kmc_count() {
    local tag="$1" input_arg="$2" db_prefix="$3" tmp_dir="$4" time_f="$5" log_f="$6"
    mkdir -p "$tmp_dir"
    /usr/bin/time -v -o "$time_f" \
        "$KMC" -k"$K" -m"$KMC_RAM_GB" -ci1 -cs4294967295 -fm -hp -t"$THREADS" \
        "$input_arg" "$db_prefix" "$tmp_dir" \
        >"$log_f" 2>&1
}

run_kmc_dump() {
    local db_prefix="$1" out_path="$2" time_f="$3" log_f="$4"
    /usr/bin/time -v -o "$time_f" \
        "$KMC_DUMP" -ci1 "$db_prefix" "$out_path" \
        >>"$log_f" 2>&1
}

run_case() {
    local scenario="$1" n="$2" rep="$3" input_arg="$4" input_bytes="$5"
    local tag="${scenario}_n${n}_r${rep}"

    echo "[bench] $tag core"
    local tuna_core_time="$RUN_DIR/${tag}.tuna_core.time"
    local tuna_core_err="$RUN_DIR/${tag}.tuna_core.stderr"
    run_tuna "${tag}.core" "$input_arg" /dev/null "$tuna_core_time" "$tuna_core_err"

    local kmc_core_time="$RUN_DIR/${tag}.kmc_count.time"
    local kmc_log="$RUN_DIR/${tag}.kmc.log"
    local kmc_db="$RUN_DIR/${tag}.kmc_db/out"
    local kmc_tmp="$RUN_DIR/${tag}.kmc_tmp"
    mkdir -p "$(dirname "$kmc_db")"
    run_kmc_count "$tag" "$input_arg" "$kmc_db" "$kmc_tmp" "$kmc_core_time" "$kmc_log"

    echo "[bench] $tag e2e"
    local tuna_e2e_time="$RUN_DIR/${tag}.tuna_e2e.time"
    local tuna_e2e_err="$RUN_DIR/${tag}.tuna_e2e.stderr"
    local tuna_out="$RUN_DIR/${tag}.tuna.tsv"
    run_tuna "${tag}.e2e" "$input_arg" "$tuna_out" "$tuna_e2e_time" "$tuna_e2e_err"

    local dump_time="$RUN_DIR/${tag}.kmc_dump.time"
    local kmc_out="$RUN_DIR/${tag}.kmc.tsv"
    run_kmc_dump "$kmc_db" "$kmc_out" "$dump_time" "$kmc_log"

    local tuna_core_wall tuna_core_rss tuna_e2e_wall tuna_e2e_rss
    tuna_core_wall=$(wall_from_time "$tuna_core_time")
    tuna_core_rss=$(rss_mb_from_time "$tuna_core_time")
    tuna_e2e_wall=$(wall_from_time "$tuna_e2e_time")
    tuna_e2e_rss=$(rss_mb_from_time "$tuna_e2e_time")

    local kmc_core_wall kmc_core_rss kmc_dump_wall kmc_dump_rss kmc_e2e_wall kmc_e2e_rss
    kmc_core_wall=$(wall_from_time "$kmc_core_time")
    kmc_core_rss=$(rss_mb_from_time "$kmc_core_time")
    kmc_dump_wall=$(wall_from_time "$dump_time")
    kmc_dump_rss=$(rss_mb_from_time "$dump_time")
    kmc_e2e_wall=$(sum_s "$kmc_core_wall" "$kmc_dump_wall")
    kmc_e2e_rss=$(max_i "$kmc_core_rss" "$kmc_dump_rss")

    local p1 p2 unique n_parts superkmers
    p1=$(parse_tuna_key "$tuna_core_err" phase1)
    p2=$(parse_tuna_key "$tuna_core_err" phase2)
    unique=$(parse_tuna_key "$tuna_core_err" unique_kmers)
    n_parts=$(parse_tuna_key "$tuna_core_err" n_parts)
    superkmers=$(parse_tuna_key "$tuna_core_err" superkmers)

    local kmc_unique
    kmc_unique=$(grep -i "unique counted k-mers\\|unique k-mers" "$kmc_log" | grep -vi "below\\|above" | awk '{print $NF}' | head -1 || true)
    [[ -n "$kmc_unique" ]] || kmc_unique=0

    echo "tuna,$scenario,$n,$rep,$tuna_core_wall,$tuna_core_rss,$tuna_e2e_wall,$tuna_e2e_rss,$p1,$p2,$unique,$n_parts,$superkmers,$input_bytes,core=/dev/null e2e=tsv" >> "$CSV"
    echo "kmc,$scenario,$n,$rep,$kmc_core_wall,$kmc_core_rss,$kmc_e2e_wall,$kmc_e2e_rss,NA,NA,$kmc_unique,NA,NA,$input_bytes,count_only_core dump_tsv_e2e" >> "$CSV"

    rm -rf "$kmc_tmp" "$(dirname "$kmc_db")"
}

generate_base

echo "[bench] results: $CSV"
echo "[bench] base: $BASE_FA"
echo "[bench] scenarios: $SCENARIOS"
echo "[bench] N values: $N_VALUES"

for scenario in $SCENARIOS; do
    for n in $N_VALUES; do
        case "$scenario" in
            repeated)
                fof=$(prepare_repeated_fof "$n")
                input_arg="@$fof"
                input_bytes=$(awk '{s += '$BASE_MB' * 1024 * 1024} END{print s}' "$fof")
                ;;
            concat)
                concat=$(prepare_concat_file "$n")
                input_arg="$concat"
                input_bytes=$(stat -c%s "$concat")
                ;;
            *)
                echo "ERROR: unknown scenario: $scenario" >&2
                exit 1
                ;;
        esac
        for rep in $(seq 1 "$REPS"); do
            run_case "$scenario" "$n" "$rep" "$input_arg" "$input_bytes"
        done
    done
done

echo "[bench] done: $CSV"
