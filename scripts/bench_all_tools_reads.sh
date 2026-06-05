#!/usr/bin/env bash
# =============================================================================
# bench_all_tools_reads.sh — k-mer counter benchmark: tuna · KMC · FastK · KFC
# Scenario : count + text dump, k=31, sequencing reads (FASTQ)
# Datasets : SRR2584863_1.fastq.gz (E. coli ~50x) and
#            SRR622461_1.fastq.gz  (Human NA12878 ~8-10x, one lane)
# Threads  : 8 (fixed)
#
# Usage: bash bench_all_tools_reads.sh <local|cluster>
#
# Notes on tool differences:
#   tuna  : single command, writes TSV directly (no dump step)
#   KMC   : count → binary, then kmc_dump → TSV  (both timed together)
#           -fq is required for FASTQ input
#   FastK : count → .ktab binary, then Tabex -A → TSV (both timed together)
#           .fastq.gz extension is recognised natively (no symlinks needed)
#   KFC   : count → .kff binary, then kfc dump → TSV (both timed together)
#           no RAM budget flag exposed
#
# Output: $WORKDIR/results.log
#   one "=== ... ===" header + /usr/bin/time -v block per run
# =============================================================================
set -uo pipefail

# =============================================================================
# Mode
# =============================================================================
MODE="${1:-}"
if [[ "$MODE" != "local" && "$MODE" != "cluster" ]]; then
    echo "Usage: $0 <local|cluster>"
    exit 1
fi

# =============================================================================
# Paths — local vs cluster
# =============================================================================
if [[ "$MODE" == "local" ]]; then
    TUNA="/home/vlevallo/documents/giulio_colab/softs/tuna/build/tuna"
    KMC="/home/vlevallo/documents/giulio_colab/softs/kmc/bin/kmc"
    KMC_DUMP="/home/vlevallo/documents/giulio_colab/softs/kmc/bin/kmc_dump"
    FASTK="/home/vlevallo/documents/giulio_colab/softs/FASTK/FastK"
    TABEX="/home/vlevallo/documents/giulio_colab/softs/FASTK/Tabex"
    KFC="/home/vlevallo/documents/giulio_colab/softs/KFC/target/release/kfc"
    DATA_DIR="$HOME/tmp/data/data_sequencing"
    WORKDIR="$HOME/tuna_bench_tmp/bench_reads_$(date +%Y%m%d_%H%M%S)"
    RAM_GB=16
else
    # On the cluster, tool paths are injected via environment variables.
    # Set them before calling this script, e.g.:
    #   export TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
    #   export KMC=/path/to/kmc  KMC_DUMP=/path/to/kmc_dump
    #   export FASTK=/path/to/FastK  TABEX=/path/to/Tabex
    #   export KFC=/path/to/kfc
    : "${TUNA:=""}"
    : "${KMC:=""}"
    : "${KMC_DUMP:=""}"
    : "${FASTK:=""}"
    : "${TABEX:=""}"
    : "${KFC:=""}"
    DATA_DIR="/WORKS/vlevallois/data/data_sequencing"
    WORKDIR="/WORKS/vlevallois/expes_tuna/bench_reads_$(date +%Y%m%d_%H%M%S)"
    RAM_GB=256
fi

ECOLI_FILE="$DATA_DIR/SRR2584863_1.fastq.gz"
HUMAN_FILE="$DATA_DIR/SRR622461_1.fastq.gz"

K=31
THREADS=8

TOOLS=(tuna kmc fastk kfc)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for bin in "$TUNA" "$KMC" "$KMC_DUMP" "$FASTK" "$TABEX" "$KFC"; do
    [[ ! -x "$bin" ]] && { echo "[error] not executable: $bin"; err=1; }
done
for f in "$ECOLI_FILE" "$HUMAN_FILE"; do
    [[ ! -f "$f" ]] && { echo "[error] file not found: $f"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

# =============================================================================
# Directory setup
# =============================================================================
for tool in tuna kmc fastk kfc; do
    mkdir -p "$WORKDIR/$tool/runs"
done

RESULTS="$WORKDIR/results.log"
{
    echo "# bench_all_tools_reads.sh — $(date)"
    echo "# MODE=$MODE"
    echo "# TUNA=$TUNA"
    echo "# KMC=$KMC  KMC_DUMP=$KMC_DUMP"
    echo "# FASTK=$FASTK  TABEX=$TABEX"
    echo "# KFC=$KFC"
    echo "# WORKDIR=$WORKDIR"
    echo "# ECOLI=$ECOLI_FILE"
    echo "# HUMAN=$HUMAN_FILE"
    echo "# K=$K  RAM=${RAM_GB}GB  THREADS=$THREADS"
} > "$RESULTS"

# =============================================================================
# Helper: write a labelled header to results.log
# =============================================================================
log_header() {
    { echo ""; echo "=== $* DATE=$(date +%Y-%m-%dT%H:%M:%S) ==="; } >> "$RESULTS"
}

# =============================================================================
# Per-tool run functions (single file)
# =============================================================================

_run_tuna() {
    local reads_file="$1" outdir="$2"
    /usr/bin/time -v \
        "$TUNA" \
            -k "$K" \
            -t "$THREADS" \
            -ram "$RAM_GB" \
            -hp \
            -w "$outdir/work" \
            "$reads_file" \
            "$outdir/out.tsv" \
        >> "$RESULTS" 2>&1
}

_run_kmc() {
    local reads_file="$1" outdir="$2"
    mkdir -p "$outdir/tmp"
    # -fq : FASTQ input (required for reads; use -fm for assembled FASTA)
    /usr/bin/time -v bash -c "
        set -e
        \"$KMC\" \
            -k${K} -ci1 -cs4294967295 -fq -m${RAM_GB} -hp -t${THREADS} \
            \"$reads_file\" \
            \"$outdir/out\" \
            \"$outdir/tmp\" \
        && \"$KMC_DUMP\" \"$outdir/out\" \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_fastk() {
    local reads_file="$1" outdir="$2"
    mkdir -p "$outdir/tmp"
    # -t1: report k-mers with count >= 1
    # .fastq.gz is recognised natively by FastK — no symlink needed
    /usr/bin/time -v bash -c "
        set -e
        \"$FASTK\" \
            -k${K} -t1 -T${THREADS} -M${RAM_GB} \
            -N\"$outdir/out\" \
            -P\"$outdir/tmp\" \
            \"$reads_file\" \
        && \"$TABEX\" -A \"$outdir/out\" > \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_kfc() {
    local reads_file="$1" outdir="$2"
    /usr/bin/time -v bash -c "
        set -e
        \"$KFC\" build \
            -k ${K} \
            -t 1 \
            -T ${THREADS} \
            -i \"$reads_file\" \
            -o \"$outdir/out.kff\" \
        && \"$KFC\" dump \
            -i \"$outdir/out.kff\" \
            -o \"$outdir/out.tsv\" \
            -t 1 \
            -T ${THREADS}
    " >> "$RESULTS" 2>&1
}

# =============================================================================
# Dispatcher: run one tool on one dataset and clean up
# =============================================================================
run_one() {
    local tool="$1" dataset="$2" reads_file="$3"
    local outdir="$WORKDIR/$tool/runs/${dataset}_t${THREADS}"
    mkdir -p "$outdir"

    log_header "TOOL=$tool DS=$dataset FILE=$(basename "$reads_file") T=$THREADS"

    case "$tool" in
        tuna)  _run_tuna  "$reads_file" "$outdir" ;;
        kmc)   _run_kmc   "$reads_file" "$outdir" ;;
        fastk) _run_fastk "$reads_file" "$outdir" ;;
        kfc)   _run_kfc   "$reads_file" "$outdir" ;;
    esac || echo "[warn] $tool failed: $dataset" | tee -a "$RESULTS"

    rm -rf "$outdir"
}

# =============================================================================
# Main loop: run all tools on each dataset
# =============================================================================
echo "[bench] Starting — $(date)"
echo "[bench] MODE=$MODE  K=$K  T=$THREADS  RAM=${RAM_GB}GB"
echo "[bench] E. coli: $(basename "$ECOLI_FILE")"
echo "[bench] Human  : $(basename "$HUMAN_FILE")"

for tool in "${TOOLS[@]}"; do
    echo "[bench]   $tool ecoli t=$THREADS — $(date +%H:%M:%S)"
    run_one "$tool" ecoli "$ECOLI_FILE"
done
for tool in "${TOOLS[@]}"; do
    echo "[bench]   $tool human t=$THREADS — $(date +%H:%M:%S)"
    run_one "$tool" human "$HUMAN_FILE"
done

echo "" >> "$RESULTS"
echo "# Done — $(date)" >> "$RESULTS"
echo "[bench] Done — results in $RESULTS"
