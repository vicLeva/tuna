#!/usr/bin/env bash
# kmer_stats.sh - k-mer count statistics for KMC databases kept on disk.
#
#   KMC_TOOLS=/path/to/kmc_tools bash kmer_stats.sh <db-prefix> [<db-prefix> ...]
#
# where <db-prefix> is the path without the .kmc_pre/.kmc_suf suffix, as left
# behind by a bench script run with KEEP_OUTPUT set.
#
# The bench scripts already write these statistics for every input they measure
# whenever KMC_TOOLS is set (see $ROOT/<experiment>_kmer_stats.csv). This script
# is for databases you kept and want to re-analyse on their own.

set -uo pipefail
export LC_ALL=C LANG=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/bench_common.sh"          # for kmer_stats_line

[[ -z "$KMC_TOOLS" || ! -x "$KMC_TOOLS" ]] && { echo "[error] set KMC_TOOLS to the kmc_tools binary" >&2; exit 1; }
(( $# )) || { echo "usage: KMC_TOOLS=... $0 <db-prefix> [...]" >&2; exit 1; }

echo "dataset,key,label,distinct_kmers,total_kmers,singletons,pct_singletons,mean_count,median_count,p99_count,max_count,cov_peak,top1pct_mass_pct,gini"
for db in "$@"; do
    db="${db%.kmc_pre}"; db="${db%.kmc_suf}"
    [[ -s "${db}.kmc_suf" ]] || { echo "[warn] no KMC db at $db" >&2; continue; }
    line=$(kmer_stats_line "$(basename "$db"),," "$db" "${db}.hist")
    [[ -n "$line" ]] && echo "$line" || echo "[warn] histogram failed for $db" >&2
done
