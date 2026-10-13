#!/usr/bin/env bash
# run_paper_suite.sh - every experiment the paper reports, in one invocation.
#
# Each bench_*.sh script measures one configuration. The paper's big-data table
# reports four of them (two k values against two memory budgets), which is why
# running the scripts once left that table incomplete. This driver runs each
# experiment over the configurations the paper actually needs, in increasing
# order of cost, and records what it ran.
#
# The steps run fastest first, so that a run stopped partway has finished as
# many experiments as possible. Measured tuna-only wall time of each, from the
# previous set of numbers:
#
#   all_datasets    0.2 h
#   scaling         0.5 h
#   thread_sweep    0.7 h
#   coverage_sweep  1.6 h   (plus building the replicated read files)
#   big_data        4.0 h   (four configurations, about an hour each)
#   msweep         11.0 h
#   nsweep         47.0 h   (two passes over 18 partition counts)
#
# Everything is resumable: each script skips rows already in its CSV, so an
# interrupted suite can be restarted and picks up where it stopped. Stopping
# after a day therefore leaves everything up to msweep complete.
#
# Required:
#   M             minimizer length. Not defaulted on purpose: it is the one
#                 setting that has to be identical across every experiment,
#                 and a stale default is what previously left the coverage
#                 sweep on a different m than the rest of the paper.
#   TUNA          tuna binary. A -DFIXED_K build only accepts its own k, so
#                 when the big-data matrix spans two k values, export
#                 TUNA_K27 and TUNA_K31 instead (or both, plus TUNA).
#   KMC, KMC_DUMP, FASTK, TABEX   only if those tools are in ONLY_TOOLS.
#
# Useful:
#   ROOT          output base, one CSV per experiment at its root
#   ONLY_TOOLS    which tools to run. Defaults to tuna alone: KMC and FastK
#                 have not changed, so their rows are merged in from the
#                 previous set afterwards instead of being measured again.
#                 Set ONLY_TOOLS="tuna kmc fastk" to measure everything.
#   ONLY_EXPES    subset of the step names below
#   BIG_KS        k values for the big-data matrix   [default: 31 27]
#   BIG_BUDGETS   memory budgets, GB, for it        [default: 256 8]
#
# Usage:
#   M=17 TUNA=/path/tuna ONLY_TOOLS=tuna bash run_paper_suite.sh

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${ROOT:=/WORKS/vlevallois/expes_tuna/paper}"
: "${BIG_KS:=31 27}"
: "${BIG_BUDGETS:=256 8}"
: "${ONLY_TOOLS:=tuna}"
: "${ONLY_EXPES:=}"

# The configuration the paper's figures read from the flat CSVs, and the one
# big-data configuration copied to $ROOT/big_data.csv for make_figures.R.
: "${K_MAIN:=31}"
: "${RAM_MAIN:=256}"

# In execution order, cheapest first.
STEPS="all_datasets scaling thread_sweep coverage_sweep big_data msweep nsweep"

step_enabled() { [[ -z "$ONLY_EXPES" || " $ONLY_EXPES " == *" $1 "* ]]; }

# Mirrors tuna_for_k in bench_common.sh, for validating before spending hours.
tuna_for_k() {
    local v="TUNA_K$1"
    [[ -n "${!v:-}" ]] && { printf '%s\n' "${!v}"; return; }
    printf '%s\n' "${TUNA:-}"
}

# --- checks ----------------------------------------------------------------
# M first and fatal: everything below reads it.
if [[ -z "${M:-}" ]]; then
    echo "[error] M is not set. It must be the same for every experiment:"
    echo "        M=17 bash run_paper_suite.sh"
    exit 1
fi

err=0

# A typo in ONLY_EXPES would otherwise run nothing at all, silently.
for want in $ONLY_EXPES; do
    [[ " $STEPS " == *" $want "* ]] || {
        echo "[error] ONLY_EXPES: unknown step '$want'"
        echo "        known steps: $STEPS"
        err=1
    }
done

wants_tuna=1
[[ -n "$ONLY_TOOLS" && " $ONLY_TOOLS " != *" tuna "* ]] && wants_tuna=0

# Which k values this run needs a tuna binary for.
needed_ks="$K_MAIN"
if step_enabled big_data; then
    for k in $BIG_KS; do
        [[ " $needed_ks " == *" $k "* ]] || needed_ks="$needed_ks $k"
    done
fi

if (( wants_tuna )); then
    for k in $needed_ks; do
        tb=$(tuna_for_k "$k")
        if [[ -z "$tb" ]]; then
            echo "[error] no tuna binary for k=$k: set TUNA_K$k (or TUNA)"; err=1
        elif [[ ! -x "$tb" ]]; then
            echo "[error] not executable: $tb (k=$k)"; err=1
        else
            # A -DFIXED_K (or -DFIXED_M) build exits on any other k or m, and
            # a default build now carries only auto_m(k), so it rejects any
            # other -m unless it was made with -DTUNA_ALL_M=ON. Prove this
            # binary accepts this (k, M) on a few hundred bases now, rather
            # than discovering it after the cheap experiments have run.
            probe=$(mktemp -d)
            { echo ">probe"
              for _ in $(seq 1 40); do printf 'ACGTTGCAGGTCATAC'; done
              echo
            } > "$probe/p.fa"
            if ! "$tb" -k "$k" -m "$M" -t 1 -co -hp -w "$probe/w" \
                    "$probe/p.fa" /dev/null >/dev/null 2>"$probe/err"; then
                echo "[error] $tb cannot count k=$k m=$M:"
                sed 's/^/        /' "$probe/err" | head -3
                err=1
            fi
            rm -rf "$probe"
        fi
    done
fi

(( err )) && exit 1

mkdir -p "$ROOT" "$ROOT/run_logs"

# --- provenance ------------------------------------------------------------
INFO="$ROOT/RUN_INFO.txt"
{
    echo "run_paper_suite.sh"
    echo "started      : $(date -Is)"
    echo "host         : $(hostname)"
    echo "git          : $(cd "$HERE/.." && git rev-parse --short HEAD 2>/dev/null) $(cd "$HERE/.." && git rev-parse --abbrev-ref HEAD 2>/dev/null)"
    echo "M            : $M"
    echo "K_MAIN       : $K_MAIN   RAM_MAIN: ${RAM_MAIN}GB"
    echo "big matrix   : k={$BIG_KS} x ram={$BIG_BUDGETS}GB"
    echo "ONLY_TOOLS   : ${ONLY_TOOLS:-all}"
    echo "ONLY_EXPES   : ${ONLY_EXPES:-all}"
    for k in $needed_ks; do
        tb=$(tuna_for_k "$k")
        echo "tuna k=$k     : $tb ($(stat -c %y "$tb" 2>/dev/null | cut -d. -f1))"
    done
    echo "ROOT         : $ROOT"
    echo
} >> "$INFO"
echo "[suite] provenance appended to $INFO"

run_step() {
    local name="$1"; shift
    local log="$ROOT/run_logs/${name}.log"
    echo ""
    echo "######## $name  ($(date +%H:%M:%S))"
    if "$@" 2>&1 | tee -a "$log"; then
        echo "[suite] $name done"
    else
        echo "[suite] $name FAILED (see $log); continuing with the rest"
    fi
}

# --- the experiments, cheapest first ---------------------------------------
# all_datasets, scaling, thread_sweep, nsweep and the coverage sweep are all
# single-configuration: one k, one budget, as the paper reports them.
for name in all_datasets scaling thread_sweep; do
    step_enabled "$name" || continue
    run_step "$name" env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        FORCE_M="$M" M="$M" ONLY_TOOLS="$ONLY_TOOLS" \
        bash "$HERE/bench_$name.sh"
done

# The coverage sweep does not source bench_common.sh, so it needs M directly.
if step_enabled coverage_sweep; then
    run_step coverage_sweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        M="$M" ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_coverage_sweep.sh"
fi

# The big-data matrix: one sub-directory per configuration, because the CSV
# name is fixed per experiment. The main configuration is also copied to
# $ROOT/big_data.csv, which is where make_figures.R looks.
if step_enabled big_data; then
    for k in $BIG_KS; do
        tb=$(tuna_for_k "$k")
        if (( wants_tuna )) && [[ ! -x "${tb:-}" ]]; then
            echo "[suite] big_data k=$k skipped: no tuna binary (set TUNA_K$k)"
            continue
        fi
        for ram in $BIG_BUDGETS; do
            sub="$ROOT/big_data_k${k}_ram${ram}"
            mkdir -p "$sub"
            run_step "big_data_k${k}_ram${ram}" env ROOT="$sub" K="$k" \
                RAM_GB="$ram" FORCE_M="$M" M="$M" TUNA="${tb:-}" \
                ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_big_data.sh"
        done
    done
    main="$ROOT/big_data_k${K_MAIN}_ram${RAM_MAIN}/big_data.csv"
    if [[ -f "$main" ]]; then
        cp -f "$main" "$ROOT/big_data.csv"
        echo "[suite] k=$K_MAIN ram=${RAM_MAIN}GB copied to \$ROOT/big_data.csv for the figures"
    fi
fi

# Tuning rather than comparison, and 11 hours of it.
if step_enabled msweep; then
    run_step msweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        bash "$HERE/bench_msweep.sh"
fi

# Last by a wide margin: two passes over 18 partition counts on the two whole
# read collections, about as long as everything above it put together, twice.
if step_enabled nsweep; then
    run_step nsweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        FORCE_M="$M" M="$M" ONLY_TOOLS="$ONLY_TOOLS" \
        bash "$HERE/bench_nsweep.sh"
fi

echo ""
echo "[suite] finished $(date -Is)"
echo "[suite] results in $ROOT"
ls -1 "$ROOT"/*.csv 2>/dev/null | sed 's/^/  /'
ls -d "$ROOT"/big_data_k*_ram* 2>/dev/null | sed 's/^/  /'
echo "finished     : $(date -Is)" >> "$INFO"
