#!/usr/bin/env bash
# scripts/bench_all_datasets_v2.sh — Rob's tuna, per-file sweep across all datasets
#
# Mirrors bench_all_datasets.sh but tuna-only (Rob's fork).
# No KMC, no correctness check, no alternating order.
# Adds G. gallus reads and H. sapiens reads (human3) as new datasets.
#
# Per-file experiment: for each dataset, run tuna on each individual file.
# Measures single-file throughput; compare with bench_all_datasets.sh output.
#
# m per dataset:
#   human assemblies : 23  (highly repetitive assembled genome)
#   everything else  : 21
#
# Output: $RESULTS_DIR/bench_datasets.csv
#   dataset,file_idx,filename,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers
# Aux files in $RESULTS_DIR/aux_datasets/

set -uo pipefail
export LC_ALL=C LANG=C

# =============================================================================
# CONFIGURE
# =============================================================================
: "${TUNA:=""}"

THREADS=${THREADS:-8}
K=31
RAM_GB=256

RESULTS_DIR="/WORKS/vlevallois/expes_tuna/rob_v2"
WORKDIR="/WORKS/vlevallois/expes_tuna"

# Dataset registry: "name:fof_path:max_files:m"
DATASETS=(
    "ecoli:/WORKS/vlevallois/data/dataset_genome_ecoli/fof.list:100:21"
    "human:/WORKS/vlevallois/data/dataset_genome_human/fof.list:10:23"
    "salmonella:/WORKS/vlevallois/data/dataset_pangenome_salmonella/fof.list:100:21"
    "gut:/WORKS/vlevallois/data/dataset_metagenome_gut/fof.list:100:21"
    "tara:/WORKS/vlevallois/data/dataset_metagenome_tara/fof.list:10:21"
    "gallus:/WORKS/vlevallois/data/dataset_reads_gallus/fof.list:999:21"
    "human3:/WORKS/vlevallois/data/dataset_reads_human3/fof.list:999:21"
)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
[[ -z "$TUNA" ]] && { echo "[error] TUNA is not set"; err=1; }
[[ -n "$TUNA" && ! -x "$TUNA" ]] && { echo "[error] not executable: $TUNA"; err=1; }
[[ "$err" -eq 1 ]] && exit 1

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$RESULTS_DIR" "$RESULTS_DIR/aux_datasets"
CSV="$RESULTS_DIR/bench_datasets.csv"
echo "dataset,file_idx,filename,m,wall_s,rss_mb,phase1_s,phase2_s,n_parts,unique_kmers" > "$CSV"

echo "[bench] Rob's tuna — per-file dataset sweep"
echo "[bench] TUNA=$TUNA"
echo "[bench] k=$K  threads=$THREADS  ram=${RAM_GB}GB"
echo "[bench] Results: $CSV"

# =============================================================================
# Helpers
# =============================================================================
wall_from_file() {
    local t; t=$(grep "Elapsed (wall clock)" "$1" | awk '{print $NF}')
    awk -F: '{if(NF==3) printf "%.3f",$1*3600+$2*60+$3; else printf "%.3f",$1*60+$2}' <<< "$t"
}
rss_mb() {
    local kb; kb=$(grep "Maximum resident" "$1" | awk '{print $NF}')
    awk "BEGIN{printf \"%.0f\", $kb/1024}"
}
se_val() { grep "^${1}:" "$2" 2>/dev/null | awk -F: '{v=$2; gsub(/^ +| +s$/,"",v); print v}' || echo na; }

# =============================================================================
# Main loop
# =============================================================================
echo ""
echo "[bench] Starting — $(date)"

for ENTRY in "${DATASETS[@]}"; do
    IFS=: read -r DS_NAME FOF DS_MAX M <<< "$ENTRY"

    if [[ ! -f "$FOF" ]]; then
        echo "  [SKIP] $DS_NAME: fof not found: $FOF"
        continue
    fi

    TOTAL=$(wc -l < "$FOF")
    N=$(( TOTAL < DS_MAX ? TOTAL : DS_MAX ))
    mapfile -t FILES < <(head -n "$N" "$FOF")

    echo ""
    echo "──── $DS_NAME  ($N / $TOTAL files, m=$M) ────────────────────"

    for IDX in "${!FILES[@]}"; do
        FILE="${FILES[$IDX]}"
        FNAME=$(basename "$FILE")
        FIDX=$(( IDX + 1 ))
        TAG="${DS_NAME}_f$(printf '%04d' "$FIDX")"

        local_tf="$RESULTS_DIR/aux_datasets/${TAG}.time"
        local_se="$RESULTS_DIR/aux_datasets/${TAG}.stderr"
        out="$WORKDIR/tuna_ds_${TAG}.tsv"
        work="$WORKDIR/tuna_ds_work_${TAG}"
        mkdir -p "$work"

        printf "  [%d/%d] %s\n" "$FIDX" "$N" "$FNAME"

        /usr/bin/time -v -o "$local_tf" \
            "$TUNA" -k "$K" -m "$M" -t "$THREADS" -ram "$RAM_GB" -hp \
            -w "$work/" "$FILE" "$out" \
            > /dev/null 2>"$local_se" \
        || { echo "    [FAIL]"; rm -rf "$work" "$out"; continue; }
        rm -rf "$work" "$out"

        wall=$(wall_from_file "$local_tf")
        rss=$(rss_mb "$local_tf")
        p1=$(se_val "phase1" "$local_se")
        p2=$(se_val "phase2" "$local_se")
        n_parts=$(se_val "n_parts" "$local_se")
        unique=$(se_val "unique_kmers" "$local_se")

        printf "    wall=%ss  p1=%ss p2=%ss  RSS=%s MB  n=%s  unique=%s\n" \
            "$wall" "$p1" "$p2" "$rss" "$n_parts" "$unique"
        echo "$DS_NAME,$FIDX,$FNAME,$M,$wall,$rss,$p1,$p2,$n_parts,$unique" >> "$CSV"
    done
done

echo ""
echo "[bench] Done — $(date)"
echo "[bench] Results: $CSV  ($(( $(wc -l < "$CSV") - 1 )) rows)"
