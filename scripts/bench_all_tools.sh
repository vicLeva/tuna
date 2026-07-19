#!/usr/bin/env bash
# =============================================================================
# bench_all_tools.sh — k-mer counter benchmark: tuna · KMC · FastK · KFC
# Scenario : count + text dump, k=31, m=21, one genome per species
# Datasets : first genome from ECOLI_FOF, first genome from HUMAN_FOF
# Threads  : 1, 4, 8, 16   —   RAM budget: 256 GB
#
# Full pipeline for every tool — no count-only/skip-output shortcuts; each
# tool writes its real output to disk, timed end to end. No phase1/phase2
# breakdown here (that's reserved for the tuna-vs-KMC scripts specifically).
#
# Notes on tool differences:
#   tuna  : single command, writes TSV directly (no dump step)
#   KMC   : count → binary, then kmc_dump → TSV  (both timed together)
#           -fm is required for assembled FASTA (multi-sequence per file)
#   FastK : count → .ktab binary, then Tabex -A → TSV (both timed together)
#           .fna extension is not recognised; symlinks to .fa are created once
#           at setup (not timed)
#   KFC   : count → .kff binary, then kfc dump → TSV (both timed together)
#           no RAM budget flag exposed
#
# Output: $WORKDIR/results.log
#   one "=== ... ===" header + /usr/bin/time -v block per run
#   bring this file back locally and grep/parse for wall time and peak RSS
# =============================================================================
set -uo pipefail

# =============================================================================
# CONFIGURE — fill in before submitting to the cluster
# (values already set in the environment take precedence)
# =============================================================================
: "${TUNA:=""}"       # /path/to/tuna
: "${KMC:=""}"        # /path/to/kmc
: "${KMC_DUMP:=""}"   # /path/to/kmc_dump
: "${FASTK:=""}"      # /path/to/FastK
: "${TABEX:=""}"      # /path/to/Tabex
: "${KFC:=""}"        # /path/to/kfc

EXPES_ROOT="/WORKS/vlevallois/expes_tuna"
WORKDIR="$EXPES_ROOT/bench_all_tools_$(date +%Y%m%d_%H%M%S)"
ECOLI_FOF="/WORKS/vlevallois/data/dataset_genome_ecoli/fof_1000.list"    # absolute path, one genome path per line, 1000 E. coli genomes
HUMAN_FOF="/WORKS/vlevallois/data/dataset_genome_human/fof_10.list"      # absolute path, one genome path per line, 10 human genomes

K=31
M=21     # fixed everywhere (standing convention going forward)
RAM_GB=256
THREADS_LIST=(1 2 4 6 8 10 16 22 28 32)

TOOLS=(tuna kmc fastk kfc)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for var in TUNA KMC KMC_DUMP FASTK TABEX KFC WORKDIR ECOLI_FOF HUMAN_FOF; do
    [[ -z "${!var}" ]] && { echo "[error] $var is not set"; err=1; }
done
for var in TUNA KMC KMC_DUMP FASTK TABEX KFC; do
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] ${!var} is not executable"; err=1; }
done
for var in ECOLI_FOF HUMAN_FOF; do
    [[ -n "${!var}" && ! -f "${!var}" ]] && { echo "[error] ${!var}: file not found"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

ECOLI_GENOME=$(head -1 "$ECOLI_FOF")
HUMAN_GENOME=$(head -1 "$HUMAN_FOF")
[[ ! -f "$ECOLI_GENOME" ]] && { echo "[error] E. coli genome not found: $ECOLI_GENOME"; exit 1; }
[[ ! -f "$HUMAN_GENOME" ]] && { echo "[error] Human genome not found: $HUMAN_GENOME"; exit 1; }

# =============================================================================
# Directory setup
# =============================================================================
mkdir -p "$WORKDIR"/fastk_links
for tool in tuna kmc fastk kfc; do
    mkdir -p "$WORKDIR/$tool/runs"
done

RESULTS="$WORKDIR/results.log"
{
    echo "# bench_all_tools.sh — $(date)"
    echo "# TUNA=$TUNA"
    echo "# KMC=$KMC  KMC_DUMP=$KMC_DUMP"
    echo "# FASTK=$FASTK  TABEX=$TABEX"
    echo "# KFC=$KFC"
    echo "# WORKDIR=$WORKDIR"
    echo "# ECOLI=$ECOLI_GENOME"
    echo "# HUMAN=$HUMAN_GENOME"
    echo "# K=$K  M=$M  RAM=${RAM_GB}GB  THREADS=${THREADS_LIST[*]}"
} > "$RESULTS"

# =============================================================================
# FastK: create .fa / .fa.gz symlinks once for all input files
# FastK determines file type from extension; .fna is not recognised.
# =============================================================================
fastk_link_path() {
    # Returns the path of the FastK-compatible symlink for a given input file.
    local f="$1"
    local base
    base=$(basename "$f")
    case "$base" in
        *.fna.gz) echo "$WORKDIR/fastk_links/${base%.fna.gz}.fa.gz" ;;
        *.fna)    echo "$WORKDIR/fastk_links/${base%.fna}.fa"       ;;
        *)        echo "$WORKDIR/fastk_links/$base"                  ;;
    esac
}

echo "[setup] Creating FastK-compatible symlinks..."
for f in "$ECOLI_GENOME" "$HUMAN_GENOME"; do
    ln -sf "$f" "$(fastk_link_path "$f")"
done

# =============================================================================
# Helper: write a labelled header to results.log
# =============================================================================
log_header() {
    { echo ""; echo "=== $* DATE=$(date +%Y-%m-%dT%H:%M:%S) ==="; } >> "$RESULTS"
}

# =============================================================================
# Per-tool run functions (single genome, single thread count)
# Each function creates an outdir, runs the full pipeline under /usr/bin/time,
# appends stdout+stderr to results.log, then removes the outdir.
# =============================================================================

_run_tuna() {
    local genome_file="$1" outdir="$2" threads="$3"
    # tuna writes TSV directly — one command, no separate dump step
    /usr/bin/time -v \
        "$TUNA" \
            -k "$K" \
            -m "$M" \
            -t "$threads" \
            -ram "$RAM_GB" \
            -hp \
            -w "$outdir/work" \
            "$genome_file" \
            "$outdir/out.tsv" \
        >> "$RESULTS" 2>&1
}

_run_kmc() {
    local genome_file="$1" outdir="$2" threads="$3"
    mkdir -p "$outdir/tmp"
    # -fm : multi-FASTA (required for assembled genomes with multiple contigs)
    # -ci1: include k-mers with count >= 1
    # -hp : hide percentage progress
    /usr/bin/time -v bash -c "
        set -e
        \"$KMC\" \
            -k${K} -ci1 -cs4294967295 -fm -m${RAM_GB} -hp -t${threads} \
            \"$genome_file\" \
            \"$outdir/out\" \
            \"$outdir/tmp\" \
        && \"$KMC_DUMP\" \"$outdir/out\" \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_fastk() {
    local genome_file="$1" outdir="$2" threads="$3"
    local link
    link=$(fastk_link_path "$genome_file")
    mkdir -p "$outdir/tmp"
    # -t1: report k-mers with count >= 1
    # -N : output prefix
    # -P : temp directory for sorting
    /usr/bin/time -v bash -c "
        set -e
        \"$FASTK\" \
            -k${K} -t1 -T${threads} -M${RAM_GB} \
            -N\"$outdir/out\" \
            -P\"$outdir/tmp\" \
            \"$link\" \
        && \"$TABEX\" -A \"$outdir/out\" > \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_kfc() {
    local genome_file="$1" outdir="$2" threads="$3"
    # -t 1 (build): superkmer solidity threshold = 1 (keep singletons)
    # -t 1 (dump) : output k-mers with count >= 1
    # KFC has no RAM budget flag
    /usr/bin/time -v bash -c "
        set -e
        \"$KFC\" build \
            -k ${K} \
            -t 1 \
            -T ${threads} \
            -i \"$genome_file\" \
            -o \"$outdir/out.kff\" \
        && \"$KFC\" dump \
            -i \"$outdir/out.kff\" \
            -o \"$outdir/out.tsv\" \
            -t 1 \
            -T ${threads}
    " >> "$RESULTS" 2>&1
}

# =============================================================================
# Dispatcher: run one tool on one genome and clean up
# =============================================================================
run_one() {
    local tool="$1" dataset="$2" genome_file="$3" threads="$4"
    local outdir="$WORKDIR/$tool/runs/${dataset}_t${threads}"
    mkdir -p "$outdir"

    log_header "TOOL=$tool DS=$dataset FILE=$(basename "$genome_file") T=$threads"

    case "$tool" in
        tuna)  _run_tuna  "$genome_file" "$outdir" "$threads" ;;
        kmc)   _run_kmc   "$genome_file" "$outdir" "$threads" ;;
        fastk) _run_fastk "$genome_file" "$outdir" "$threads" ;;
        kfc)   _run_kfc   "$genome_file" "$outdir" "$threads" ;;
    esac || echo "[warn] $tool failed: $dataset t=$threads" | tee -a "$RESULTS"

    rm -rf "$outdir"
}

# =============================================================================
# Main loop: for each thread count, run all tools on each genome
# =============================================================================
echo "[bench] Starting — $(date)"
echo "[bench] E. coli: $(basename "$ECOLI_GENOME")"
echo "[bench] Human  : $(basename "$HUMAN_GENOME")"

for threads in "${THREADS_LIST[@]}"; do
    echo "[bench] ── threads=$threads ──"
    for tool in "${TOOLS[@]}"; do
        echo "[bench]   $tool ecoli t=$threads — $(date +%H:%M:%S)"
        run_one "$tool" ecoli "$ECOLI_GENOME" "$threads"
    done
    for tool in "${TOOLS[@]}"; do
        echo "[bench]   $tool human t=$threads — $(date +%H:%M:%S)"
        run_one "$tool" human "$HUMAN_GENOME" "$threads"
    done
done

echo "" >> "$RESULTS"
echo "# Done — $(date)" >> "$RESULTS"
echo "[bench] Done — results in $RESULTS"
