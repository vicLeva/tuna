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
# many experiments as possible. Measured *tuna-only* wall time of each, from
# the previous set of numbers. With KMC and FastK as well, expect roughly
# three times these, and more where FastK runs to its timeout:
#
#   all_datasets    0.4 h   (per file, then each collection as a whole)
#   thread_sweep    0.7 h
#   coverage_sweep  1.6 h   (plus building the replicated read files)
#   big_data        4.0 h   (four configurations, about an hour each)
#   nsweep         47.0 h   (two passes over 18 partition counts)
#
# Everything is resumable: each script skips rows already in its CSV, so an
# interrupted suite can be restarted and picks up where it stopped. Stopping
# after a day therefore leaves everything up to msweep complete.
#
# Required:
#   (m is not a setting here: tuna derives it from k, which is what a user
#   gets and so what the paper measures. A plain build is what you want --
#   no -DTUNA_ALL_M, no -DFIXED_K.)
#   TUNA          tuna binary. A -DFIXED_K build only accepts its own k, so
#                 when the big-data matrix spans two k values, export
#                 TUNA_K27 and TUNA_K31 instead (or both, plus TUNA).
#   KMC, KMC_DUMP, FASTK, TABEX   required unless ONLY_TOOLS excludes them.
#
# Useful:
#   ROOT          output base, one CSV per experiment at its root
#   ONLY_TOOLS    which tools to run. All three by default. Set it to "tuna"
#                 to refresh tuna alone and merge the competitor rows in from
#                 an earlier set afterwards.
#   ONLY_EXPES    subset of the step names below
#   BIG_KS        k values for the big-data matrix   [default: 31 27]
#   BIG_BUDGETS   memory budgets, GB, for it        [default: 256 8]
#
# Usage:
#   TUNA=/path/tuna KMC=... KMC_DUMP=... FASTK=... TABEX=... bash run_paper_suite.sh

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${ROOT:=/WORKS/vlevallois/expes_tuna/paper}"
: "${BIG_KS:=31 27}"
: "${BIG_BUDGETS:=256 8}"
: "${ONLY_TOOLS:=}"
: "${ONLY_EXPES:=}"

# The configuration the paper's figures read from the flat CSVs, and the one
# big-data configuration copied to $ROOT/big_data.csv for make_figures.R.
: "${K_MAIN:=31}"
: "${RAM_MAIN:=256}"

# In execution order, cheapest first.
#
# scaling and msweep are deliberately absent from the default order. The
# whole-collection pass of all_datasets replaced scaling in the paper, and the
# m sweep has served its purpose now that tuna picks m itself. Both scripts
# still work; ask for them by name with ONLY_EXPES.
STEPS="all_datasets thread_sweep coverage_sweep big_data nsweep scaling msweep"

step_enabled() { [[ -z "$ONLY_EXPES" || " $ONLY_EXPES " == *" $1 "* ]]; }

# Mirrors tuna_for_k in bench_common.sh, for validating before spending hours.
tuna_for_k() {
    local v="TUNA_K$1"
    [[ -n "${!v:-}" ]] && { printf '%s\n' "${!v}"; return; }
    printf '%s\n' "${TUNA:-}"
}

# --- checks ----------------------------------------------------------------
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
            # A -DFIXED_K (or -DFIXED_M) build exits on any other k. Prove
            # this binary accepts this k on a few hundred bases now, rather
            # than discovering it after the cheap experiments have run.
            probe=$(mktemp -d)
            { echo ">probe"
              for _ in $(seq 1 40); do printf 'ACGTTGCAGGTCATAC'; done
              echo
            } > "$probe/p.fa"
            if ! "$tb" -k "$k" -t 1 -co -hp -w "$probe/w" \
                    "$probe/p.fa" /dev/null >/dev/null 2>"$probe/err"; then
                echo "[error] $tb cannot count k=$k:"
                sed 's/^/        /' "$probe/err" | head -3
                err=1
            fi
            rm -rf "$probe"
        fi
    done
fi

# Same for the competitors: better to hear about a missing binary now than
# once per step, after the first experiment has already started.
tool_wanted() { [[ -z "$ONLY_TOOLS" || " $ONLY_TOOLS " == *" $1 "* ]]; }
for pair in "kmc:KMC KMC_DUMP" "fastk:FASTK TABEX"; do
    tool_wanted "${pair%%:*}" || continue
    for v in ${pair#*:}; do
        if [[ -z "${!v:-}" ]]; then
            echo "[error] $v is not set (needed by ${pair%%:*})"; err=1
        elif [[ ! -x "${!v}" ]]; then
            echo "[error] not executable: ${!v} ($v)"; err=1
        fi
    done
done

(( err )) && exit 1

mkdir -p "$ROOT" "$ROOT/run_logs"

# --- provenance ------------------------------------------------------------
INFO="$ROOT/RUN_INFO.txt"
{
    echo "run_paper_suite.sh"
    echo "started      : $(date -Is)"
    echo "host         : $(hostname)"
    echo "git          : $(cd "$HERE/.." && git rev-parse --short HEAD 2>/dev/null) $(cd "$HERE/.." && git rev-parse --abbrev-ref HEAD 2>/dev/null)"
    echo "m            : auto (tuna derives it from k)"
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
# all_datasets, thread_sweep, nsweep and the coverage sweep are all
# single-configuration: one k, one budget, as the paper reports them.
for name in all_datasets thread_sweep; do
    step_enabled "$name" || continue
    run_step "$name" env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_$name.sh"
done

# The coverage sweep does not source bench_common.sh, so it needs M directly.
if step_enabled coverage_sweep; then
    run_step coverage_sweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_coverage_sweep.sh"
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
                RAM_GB="$ram" TUNA="${tb:-}" \
                ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_big_data.sh"
        done
    done
    main="$ROOT/big_data_k${K_MAIN}_ram${RAM_MAIN}/big_data.csv"
    if [[ -f "$main" ]]; then
        cp -f "$main" "$ROOT/big_data.csv"
        echo "[suite] k=$K_MAIN ram=${RAM_MAIN}GB copied to \$ROOT/big_data.csv for the figures"
    fi
fi

# Superseded by all_datasets' whole-collection pass, so it is not in the
# default order; still runs if ONLY_EXPES asks for it.
if step_enabled scaling && [[ -n "$ONLY_EXPES" ]]; then
    run_step scaling env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_scaling.sh"
fi

# Tuning rather than comparison, 11 hours of it, and the m it was choosing is
# now chosen by tuna. Out of the default order; still runs on request.
if step_enabled msweep && [[ -n "$ONLY_EXPES" ]]; then
    run_step msweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        bash "$HERE/bench_msweep.sh"
fi

# Last by a wide margin: two passes over 18 partition counts on the two whole
# read collections, about as long as everything above it put together, twice.
if step_enabled nsweep; then
    run_step nsweep env ROOT="$ROOT" K="$K_MAIN" RAM_GB="$RAM_MAIN" \
        ONLY_TOOLS="$ONLY_TOOLS" bash "$HERE/bench_nsweep.sh"
fi

echo ""
echo "[suite] finished $(date -Is)"
echo "[suite] results in $ROOT"
ls -1 "$ROOT"/*.csv 2>/dev/null | sed 's/^/  /'
ls -d "$ROOT"/big_data_k*_ram* 2>/dev/null | sed 's/^/  /'
echo "finished     : $(date -Is)" >> "$INFO"
