#!/usr/bin/env bash
# Synthetic GRCh38-derived read-shaped Tuna vs KMC scaling benchmark.
#
# Generates one sampled FASTA read file per requested read count. Reads are
# sampled deterministically from a GRCh38-derived ACTG sequence so the observed
# k-mer repetition is genome-like rather than iid random DNA.
#
# Timing reported:
#   core: Tuna to /dev/null; KMC count only
#   e2e:  Tuna TSV to scratch; KMC count + kmc_dump TSV to scratch
#
# All generated data, outputs, logs, and temporary files stay under WORK.

set -euo pipefail

REF=${REF:-/scratch2/tmp/grch38_1019/genome.fa.gz}
WORK=${WORK:-/scratch3/tmp/tuna-optimization-codex/bench/synthetic_grch38_reads}
TMPDIR=${TMPDIR:-/scratch3/tmp/tuna-optimization-codex/tmp}

TUNA=${TUNA:-/scratch3/tmp/tuna-optimization-codex/build/tuna-k31-m21/tuna}
KMC=${KMC:-/scratch3/rob/tuna-optimization/KMC/bin/kmc}
KMC_DUMP=${KMC_DUMP:-/scratch3/rob/tuna-optimization/KMC/bin/kmc_dump}

K=${K:-31}
M=${M:-21}
THREADS=${THREADS:-8}
PARTS=${PARTS:-}
KMC_RAM_GB=${KMC_RAM_GB:-64}
GENOME_MB=${GENOME_MB:-256}
READ_LEN=${READ_LEN:-150}
READ_VALUES=${READ_VALUES:-"100000 250000 500000 1000000 2000000"}
SEED=${SEED:-1}
REPS=${REPS:-1}

mkdir -p "$WORK" "$TMPDIR"
export TMPDIR

DATA_DIR="$WORK/data"
RUN_DIR="$WORK/runs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$DATA_DIR" "$RUN_DIR"

GENOME_SEQ="$DATA_DIR/grch38_${GENOME_MB}mb.actg"
CSV="$RUN_DIR/bench_synthetic_grch38_reads.csv"

require_exe() {
    local p="$1"
    [[ -x "$p" ]] || { echo "ERROR: not executable: $p" >&2; exit 1; }
}

require_exe "$TUNA"
require_exe "$KMC"
require_exe "$KMC_DUMP"
[[ -f "$REF" ]] || { echo "ERROR: missing reference: $REF" >&2; exit 1; }

echo "tool,reads,rep,read_len,core_wall_s,core_rss_mb,e2e_wall_s,e2e_rss_mb,phase1_s,phase2_s,unique_kmers,n_parts,superkmers,input_bytes,notes" > "$CSV"

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

generate_genome_seq() {
    if [[ -s "$GENOME_SEQ" ]]; then
        return
    fi
    echo "[setup] generating $GENOME_SEQ from $REF (${GENOME_MB} MiB ACTG)"
    python3 - "$REF" "$GENOME_SEQ" "$GENOME_MB" <<'PY'
import gzip
import sys

ref, out, mb = sys.argv[1], sys.argv[2], int(sys.argv[3])
target = mb * 1024 * 1024
written = 0

with gzip.open(ref, "rt", encoding="ascii", errors="ignore") as inp, open(out, "w") as out_f:
    for line in inp:
        if written >= target:
            break
        if line.startswith(">"):
            continue
        seq = ''.join(ch for ch in line.upper() if ch in "ACGT")
        if not seq:
            continue
        take = min(len(seq), target - written)
        out_f.write(seq[:take])
        written += take

if written < target:
    raise SystemExit(f"reference ended after {written} bases, wanted {target}")
PY
}

prepare_reads_file() {
    local reads="$1"
    local out="$DATA_DIR/grch38_reads_len${READ_LEN}_n${reads}.fa"
    if [[ ! -s "$out" ]]; then
        echo "[setup] creating read input $out" >&2
        python3 - "$GENOME_SEQ" "$out" "$reads" "$READ_LEN" "$SEED" <<'PY'
import random
import sys

genome_path, out_path = sys.argv[1], sys.argv[2]
n_reads, read_len, seed = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
seq = open(genome_path, "r", encoding="ascii").read().strip()
if len(seq) < read_len:
    raise SystemExit("source sequence is shorter than read length")
rng = random.Random(seed)

with open(out_path, "w", encoding="ascii") as out:
    max_start = len(seq) - read_len
    for i in range(n_reads):
        start = rng.randint(0, max_start)
        out.write(f">read_{i}_{start}\n")
        out.write(seq[start:start + read_len] + "\n")
PY
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
    local input_arg="$1" db_prefix="$2" tmp_dir="$3" time_f="$4" log_f="$5"
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
    local reads="$1" rep="$2" input_arg="$3" input_bytes="$4"
    local tag="reads_n${reads}_r${rep}"

    echo "[bench] $tag core"
    local tuna_core_time="$RUN_DIR/${tag}.tuna_core.time"
    local tuna_core_err="$RUN_DIR/${tag}.tuna_core.stderr"
    run_tuna "${tag}.core" "$input_arg" /dev/null "$tuna_core_time" "$tuna_core_err"

    local kmc_core_time="$RUN_DIR/${tag}.kmc_count.time"
    local kmc_log="$RUN_DIR/${tag}.kmc.log"
    local kmc_db="$RUN_DIR/${tag}.kmc_db/out"
    local kmc_tmp="$RUN_DIR/${tag}.kmc_tmp"
    mkdir -p "$(dirname "$kmc_db")"
    run_kmc_count "$input_arg" "$kmc_db" "$kmc_tmp" "$kmc_core_time" "$kmc_log"

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

    echo "tuna,$reads,$rep,$READ_LEN,$tuna_core_wall,$tuna_core_rss,$tuna_e2e_wall,$tuna_e2e_rss,$p1,$p2,$unique,$n_parts,$superkmers,$input_bytes,core=/dev/null e2e=tsv" >> "$CSV"
    echo "kmc,$reads,$rep,$READ_LEN,$kmc_core_wall,$kmc_core_rss,$kmc_e2e_wall,$kmc_e2e_rss,NA,NA,$kmc_unique,NA,NA,$input_bytes,count_only_core dump_tsv_e2e" >> "$CSV"

    rm -rf "$kmc_tmp" "$(dirname "$kmc_db")"
}

generate_genome_seq

echo "[bench] results: $CSV"
echo "[bench] genome source: $GENOME_SEQ"
echo "[bench] read length: $READ_LEN"
echo "[bench] read counts: $READ_VALUES"

for reads in $READ_VALUES; do
    input=$(prepare_reads_file "$reads")
    input_bytes=$(stat -c%s "$input")
    for rep in $(seq 1 "$REPS"); do
        run_case "$reads" "$rep" "$input" "$input_bytes"
    done
done

echo "[bench] done: $CSV"
