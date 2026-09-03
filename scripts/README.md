# scripts/ — benchmark suite for the tuna paper

These are the scripts the published numbers come from. Earlier, ad-hoc versions
live in `benchmark/legacy/` and are kept only for reference.

## Layout

| file | what it measures |
|---|---|
| `bench_common.sh` | shared machinery: tool runners, timing, guards. Sourced, not run. |
| `bench_all_datasets.sh` | one file at a time, across five collections (tuna, KMC) |
| `bench_scaling.sh` | n input files, increasing n (tuna, KMC, FastK) |
| `bench_big_data.sh` | the two large read sets, whole (tuna, KMC, FastK) |
| `bench_thread_sweep.sh` | thread sweep T=1..32, assemblies and reads (tuna, KMC, FastK) |
| `bench_coverage_sweep.sh` | coverage 1x-100x on one read set (tuna, KMC) — **already run** |
| `bench_nsweep.sh` | partition count `-n` = 2^5..2^21: timings, plus table stats via `-dbg` at a subset of `n` (tuna only) |
| `kmer_stats.sh` | count statistics from a kept KMC database |

## Running

```bash
export TUNA=/WORKS/vlevallois/softs/tuna/build/tuna
export KMC=.../kmc  KMC_DUMP=.../kmc_dump
export FASTK=.../FastK  TABEX=.../Tabex

bash bench_thread_sweep.sh       # each experiment stands alone
bash bench_scaling.sh
...
```

Results land in `$ROOT`, default `/WORKS/vlevallois/expes_tuna/expes_paper`,
one CSV per experiment.

Useful overrides: `ROOT`, `DATA_ROOT`, `K`, `M`, `THREADS`, `RAM_GB`,
`TIMEOUT_S`, `FASTK_TIMEOUT_S`, `THREAD_LIST`, and for side experiments
`FORCE_M` (one minimizer length everywhere, ignoring the per-dataset values),
`ONLY_TOOLS` (a subset of `tuna kmc fastk`), `ONLY_MODES` (a subset of
`bin ascii`), `MAX_FILES` (cap per dataset) and
`KEEP_OUTPUT` (see below).

## k-mer statistics

Set `KMC_TOOLS` and every experiment also describes each input it measures, in
`$ROOT/<experiment>_kmer_stats.csv`:

```
dataset,key,label,distinct_kmers,total_kmers,singletons,pct_singletons,
mean_count,median_count,p99_count,max_count,cov_peak,top1pct_mass_pct,gini
```

- `singletons` — k-mers seen exactly once: mostly sequencing error in reads,
  near zero in assemblies
- `cov_peak` — most common count above 1, i.e. the coverage mode
- `top1pct_mass_pct`, `gini` — how skewed the count distribution is

This costs a few seconds per point and no extra disk. The numbers come from the
count histogram of the database the KMC run just built, taken before that
database is deleted; the histogram is a few hundred kB whatever the input size,
while an ASCII dump would be ~40 bytes per distinct k-mer — hundreds of GB on
the big read sets — and a pass over it would take hours to yield the same
figures. Histograms are kept under `$ROOT/hist/<experiment>/`.

Statistics need a KMC run, so they are not produced when `ONLY_TOOLS=tuna`.

`max_count` and the quantiles are capped at 10^6 occurrences (`-cx1000000`),
which is far above anything these datasets reach but matters if the number is
ever quoted as an absolute.

### Keeping outputs

Runs delete their output as soon as it is measured. `KEEP_OUTPUT=<dir>`
preserves the **binary** results instead — the KMC database and tuna's KFF:

```bash
KEEP_OUTPUT=$ROOT/kept bash bench_big_data.sh
KMC_TOOLS=.../kmc_tools bash kmer_stats.sh $ROOT/kept/*.kmc_pre
```

ASCII output is never kept, for the size reason above. `kmer_stats.sh` re-runs
the same analysis on databases kept this way.

Every experiment is independent and **resumable**: re-running a script only
fills the gaps in its own CSV, so one experiment can be redone without
disturbing the others.

## The two output modes

Timings are only comparable within a mode.

| mode | meaning |
|---|---|
| `bin` | the tool's own binary output: tuna KFF, KMC `.kmc_pre`/`.kmc_suf`, FastK `.ktab` |
| `ascii` | a plain-text k-mer/count table |

The tools reach ASCII differently, which is why they are timed differently:

- **tuna** writes ASCII directly, in one pass, so it gets **two separate runs**
  (one KFF, one TSV). Its ASCII cost is `ascii_wall - bin_wall`.
- **KMC** writes a binary DB and `kmc_dump` converts it: one KMC run plus a
  separately timed dump, `ascii_wall = kmc + dump`.
- **FastK** writes a `.ktab` and `Tabex -A` converts it: same shape as KMC.

So `dump_s` holds the conversion step alone for KMC and FastK, and is empty for
tuna, where the equivalent quantity is the difference between its two runs.

There is deliberately no count-only mode: two measurements per tool are enough
to isolate the price of ASCII, and a third would add ~50 % machine time.

## CSV schema

```
dataset,key,label,tool,mode,wall_s,rss_mb,phase1_s,phase2_s,dump_s,
unique_kmers,total_kmers,status
```

- `key` — the experiment's independent variable: file index (`all_datasets`),
  number of input files (`scaling`, `big_data`), thread count (`thread_sweep`).
- `label` — free text: the file name in `all_datasets`, empty elsewhere.
- `phase1_s` / `phase2_s` — tuna's own phases; KMC's 1st/2nd stage; empty for
  FastK, which reports no equivalent breakdown.
- `status` — `ok`, `fail`, or `timeout`.

Empty means *not measured*, never zero.

## Failure handling

Every run is wrapped in `timeout`. FastK is known to **hang** rather than fail
on certain file-count/thread-count combinations — it spins in `Scan_All_Input`
and never returns — so hangs are capped and recorded as `timeout`, while any
other non-zero exit is `fail`. Keeping the two distinct matters: a hang and a
crash say different things about a tool.

`Tabex` crashes on a large share of FastK tables. When it does, the `bin` row
still stands and only the `ascii` row is marked failed.

Failed rows are written, not dropped. "FastK cannot do this" is a result.

## bench_coverage_sweep.sh

Kept for provenance, **not** to be re-run: it is the exact script behind
`benchmark/results/coverage_sweep/`, produced 2026-08-28. It stands apart from
the rest on purpose — it does not source `bench_common.sh`, and its default
`ROOT` still points at `expes_tuna/coverage_sweep` — because editing it would
break the correspondence between the script and those results.

It also differs in design, for reasons specific to that experiment: one full
pass per tool rather than per input, and level files built incrementally
(`mv` + append one source copy), so each level is written once and only one
level file exists at a time.

## bench_nsweep.sh

Partition-count sweep on the two large read sets (`gallus`, `human3`), tuna
only. `-n` is rounded up to a power of two by tuna itself, so
`n = 2^5, 2^6, ..., 2^21` (plus `auto`, the value tuna picks with no `-n`) is
the entire reachable search space — nothing finer exists between two
consecutive powers of two.

Unlike every other script here, this takes **one tuna binary per run** with no
tool comparison: the question is how tuna's own partition count trades off
against its own table, not tuna vs KMC vs FastK. To compare branches (say
kache-hash vs unordered_dense), run it once per build with a different `TUNA`
and `ROOT` each time and diff the CSVs afterwards.

Each `n` gets two passes, writing two different CSVs, because they answer
different questions:

- **pass 1, timing** — every `n` in `$NS`. Measures wall time, phase1/phase2,
  and peak RSS. Output is always `.kff` (bin): the sweep is about
  partitioning and counting cost, not output format. Written to
  `$ROOT/nsweep.csv`.
- **pass 2, table stats** — every `n` in `$DBG_NS` (a coarser subset of `$NS`
  by default), run with `-dbg`. This makes tuna emit `debug_table_stats.csv`
  (one row per partition: `init_sz`, `table_cap`, `n_inserted`, `n_unique`,
  `load_factor`, `n_resizes`) and `debug_resize_events.csv` (one row per
  resize), copied to `$ROOT/dbg/<dataset>_n<n>/` and aggregated into
  `$ROOT/nsweep_parts.csv` (mean/median/min/max unique k-mers per partition,
  mean init size and capacity, mean load factor, total resizes). `-dbg`'s
  bookkeeping perturbs wall time, so **pass 2's own timings are not comparable
  to pass 1's** — it exists to look inside the tables, not to race them.
  `$DBG_NS` defaults to a subset of `$NS` because the quantities it captures
  move smoothly between points, so sampling every `n` a second time would
  double an already long sweep for no signal.

Table layout is branch-specific — kache-hash quotients the key out of the
bucket address, unordered_dense and absl::flat_hash_map store it in full — so
`table_cap` means whatever that branch calls capacity/bucket_count, and
converting it to bytes per partition needs a branch-aware constant applied at
analysis time, not inside the script. A branch whose table does not track
resizes (unordered_dense's `resize_log()` is a stub) will always show
`n_resizes = 0`: that is a property of the table, not a failed measurement.

Resumable per `(dataset, n)`, independently for each pass, same as everything
else here. Default grid is 17 points x 2 datasets = 34 timing runs, of which 8
x 2 = 16 also get a stats run — narrow `NS` or `DBG_NS` for a faster first
signal; the full run will skip whatever is already in the CSVs.

## Input regimes in the thread sweep

`bench_thread_sweep.sh` covers both assemblies and reads in one run, because
they stress a counter in opposite ways and a tool can scale well on one and
badly on the other:

- **assemblies** (`ecoli`, `human`) — error-free and non-redundant, so nearly
  every k-mer is distinct and seen once or twice. Passed as a fof, so
  parallelism comes from spreading *files* across workers.
- **reads** (`ecoli_reads`, `human_reads`) — roughly half the k-mers are
  sequencing errors seen exactly once, which inflates the table with junk and
  shifts the cost into counting. Passed as *one* compressed file, so
  parallelism has to come from inside a single stream, which is a different
  code path.

The `dataset` column tells them apart. Note that for the single-file read
inputs, tuna takes its single-`.gz` producer/consumer path, whose shape depends
on the thread count — worth remembering when reading the T=1 and T=2 points.

## Notes

- **KFC** was part of earlier thread sweeps and has been dropped.
- **FastK input names**: FastK dispatches on file extension and does not accept
  `.fna`, so inputs are exposed through a symlink farm under `$ROOT/fastk_links`.
  The real data is never touched.
- **Minimizer length** is `m=21` everywhere. A whole-experiment comparison of
  m=17 against m=21 on `all_datasets` and `big_data` found the two within noise
  of each other, so the simpler uniform value was kept.
