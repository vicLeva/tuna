# Optimization Round 2

This document records the engineering analysis, changes, and correctness gates
for the `optimization-round2` branch.

## Correctness contract

A performance result is accepted only when Tuna and KMC agree exactly on:

- total canonical k-mers, including multiplicity;
- distinct canonical k-mers after `-ci`/`-cx` filtering.

Small inputs also compare every sorted `(k-mer, count)` row. Large timing runs
use `benchmark/compare_kmc_counts.sh`; Tuna uses `-co` and KMC uses `-w`, so
neither tool pays text-output or database-serialization costs. Real-output
runs use `benchmark/compare_kmc_binary.sh`: Tuna writes KFF and KMC writes its
native binary database.

## Engineering overview

Both counters use two broad phases, but their counting engines remain
fundamentally different:

| Phase | Tuna | KMC 3 |
|---|---|---|
| 1 | Parse FASTA/FASTQ, form minimizer superkmers, and route packed records into partitions | Parse input, form signature-based superkmers, and route compact records into bins |
| 2 | Aggregate repeated packed superkmers when profitable, then increment canonical k-mers in a partition-local streaming hash table | Radix-sort compact k-mer records in each bin and compact equal runs |
| Output | Iterate hash tables and write TSV/KFF; `-co` skips serialization | Store a KMC database; `-w` skips database output |

The initial large-data profile showed two separate problems. Phase 1 repeatedly
packed overlapping superkmers, allocated one string per parser chunk, and did
gzip decompression and parsing serially for a single file. Phase 2 routed table
buckets by the superkmer minimizer. That creates excellent locality within a
superkmer, but on large and repetitive inputs it also concentrates unrelated
k-mers into too few buckets and sends billions of entries to the overflow
table.

## Implemented changes

1. **Separate partitioning from table routing.** Phase 1 now uses an `m=17`
   minimizer to make compact superkmers without sacrificing partition entropy.
   Phase 2 routes every canonical
   k-mer by Tuna's rolling k-mer hash, independently of that minimizer. This
   retains Tuna's hash-table counting method while avoiding minimizer-induced
   overflow skew.
2. **Pack each parsed sequence once.** Superkmer writers copy base-aligned
   slices from a pre-packed sequence instead of repacking every overlapping
   record.
3. **Use compact partition records.** A record contains only its length and
   packed bases. The obsolete stored minimizer coordinate was removed; the
   common `k=31, m=17` record header is one byte.
4. **Pipeline compressed input.** A dedicated gzip worker inflates into a
   reusable buffer ring while the FASTX parser fills recycled contiguous
   batches and consumer workers perform minimizer extraction. Multiple gzip
   inputs can use multiple producer pipelines without exceeding `-t`.
5. **Recycle parser batches.** ACTG chunks are stored in a reusable arena plus
   offsets, eliminating one heap allocation per chunk.
6. **Initialize the next counting window once.** Two rolling k-mer windows
   alternate, allowing the next bucket to be prefetched while the current
   superkmer is counted.
7. **Use a serial overflow table where ownership is serial.** A disk-mode
   partition is counted by one worker, so its overflow table no longer pays
   atomic claim and resize-coordination costs.
8. **Make count-only mode a correctness gate.** `-co` skips serialization but
   still reports filtered distinct and total counts, and the benchmark helper
   compares them with KMC before accepting a timing.
9. **Aggregate repeated superkmers exactly.** Phase 2 stores stable
   `string_view` keys into the mmap or in-memory partition buffer and replays
   each unique record with its multiplicity. `-dedup auto` samples redundancy
   per partition and falls back to direct counting before any table insert if
   aggregation is not profitable or exceeds its memory allowance.
10. **Normalize packed record tails.** Unused bits in the final packed byte are
    cleared so identical superkmers have identical keys regardless of adjacent
    source bases.
11. **Stabilize table sizing.** Calibrated tables use a 64K minimum initial
    size. This removes a first-completed-partition race that changed chicken
    overflow traffic from about 55 thousand to 5.1 million insertions.
12. **Batch and compact KFF output.** Counting workers build complete KFF raw
    records and submit megabyte batches, replacing two library calls and two
    copies per k-mer. Each batch uses the smallest lossless 1–4 byte count
    width by opening a standards-compliant KFF raw section when the width
    changes.
13. **Avoid duplicate rolling-hash state.** Counting already retains the
    outgoing base in its directed k-mer window. A known-outgoing-base advance
    path now updates ntHash directly instead of also maintaining a 31-byte
    circular base queue.
14. **Prefetch aggregated records.** Exact superkmer replay alternates two
    initialized windows so the next unique record's table bucket is prefetched
    while the current record is inserted.
15. **Remove KFF per-record vector growth.** KFF workers allocate one fixed
    batch and reuse it by index. Compile-time 1-, 2-, and 3-byte compaction
    loops replace runtime-sized moves in the billion-record output path.
16. **Overlap temporary-file cleanup.** Disk workers unlink each consumed
    partition after closing its mmap, distributing cleanup across Phase 2
    instead of serially unlinking 32,768 files after all workers join.
17. **Treat empty partitions as valid.** The mmap reader accepts a zero-byte
    partition as an empty stream. This is required for valid sparse minimizer
    distributions and was exposed by the minimizer sweep.

The shorter partition minimizer and rolling-hash route are a hybrid, not a KMC
counting clone: KMC sorts records, while Tuna still performs canonical rolling
hash insertion and increment semantics.

## Validation

All measurements below use canonical `k=31`, eight requested threads, an 8 GB
RAM limit, and count-only output.

### Exact-output sample

A 100,000-read sample contains 6,979,485 total and 6,387,572 distinct k-mers.
Tuna and KMC match all 6,387,572 sorted `(k-mer, count)` rows byte for byte.

### Chicken FASTQ

The first large-file check is a 1.6 GB gzipped chicken FASTQ:

| Tool/version | Phase 1 | Phase 2 | Wall | Total | Distinct |
|---|---:|---:|---:|---:|---:|
| KMC | 7.31 s | 5.13 s | 12.47 s | 1,691,878,691 | 832,852,775 |
| Tuna branch point | 10.43 s | 7.01 s | 17.79 s | 1,691,878,691 | 832,852,775 |
| Tuna optimized | 5.85 s | 4.77 s | 10.81 s | 1,691,878,691 | 832,852,775 |

The optimized current source is about 13 percent faster than KMC wall-clock on
this run and matches both count invariants exactly.

### Repeated-superkmer cost

On a larger chicken read pair with 4,425,430,316 total and 1,290,334,589
distinct k-mers, direct Phase 2 had a 12.13 s two-run mean. The final adaptive
path measured 6.89, 6.95, and 7.18 s while preserving both counts.

An instrumented run aggregated 545,908,037 records into 312,092,826 unique
records (mean multiplicity 1.749). Across eight workers, exact map construction
used 6.97 CPU-seconds and unique-record replay/counting used 47.62 CPU-seconds.
The dedup front-end is therefore about 12.8 percent of those two CPU stages,
while avoiding 42.8 percent of record replays.

### Human paired FASTQ

The human-scale gate uses the 29 GB compressed ERR3239276 read pair:

| Tool/version | Phase 1 | Phase 2 | Wall | Peak RSS | Total | Distinct |
|---|---:|---:|---:|---:|---:|---:|
| KMC | 85.81 s | 323.79 s | 415.89 s | 7,925,404 KB | 95,095,613,344 | 7,589,178,026 |
| Tuna before dev integration | 133.89 s | 210.33 s | 352.45 s | 2,547,252 KB | 95,095,613,344 | 7,589,178,026 |
| Tuna hybrid | 158.56 s | 88.03 s | 254.67 s | 2,771,052 KB | 95,095,613,344 | 7,589,178,026 |

The hybrid is about 39 percent faster wall-clock than KMC while using about 65
percent less peak RSS. Its counting phase is 73 percent faster than KMC and 58
percent faster than the pre-integration Tuna run. Phase 1 remains slower than
KMC and varied materially between the two Tuna runs, so it remains the primary
optimization target.

Timings are single-run development measurements on the same host and should be
replaced by repeated medians before making release claims.

## Dev branch audit

The useful `dev` idea was exact repeated-superkmer aggregation. The hybrid keeps
this branch's compact one-header record, shorter partition minimizer,
rolling-hash table routing, direct packed k-mer initialization, and faster Phase
1 pipeline, while adopting multiplicity replay through a lower-overhead
`string_view` map.

The complete `dev` Phase 1 was not ported: its recorded human Phase 1 was
149.19 s versus 133.89 s for the then-current pipeline. Its packed multi-file
scheduler and async writer are tightly coupled, so neither is a demonstrated
isolated win here. Compiling out minimizer-position tracking was evaluated
separately and rejected after an alternating A/B measured 10.55 s versus 10.21
s with tracking enabled.

## Experiments rejected

- `m=9` made even longer superkmers but reduced partition entropy enough to
  create badly imbalanced phase-2 work.
- Routing the table by shorter minimizers greatly increased overflow traffic.
- Enlarging per-partition write buffers increased memory and temporary storage
  without a repeatable wall-time gain.
- Leaving unused packed tail bits unnormalized saved about 2.9 percent in an
  isolated Phase 1 test but raised Phase 2 from 6.89-7.18 s to 7.34-7.37 s,
  losing end to end.
- A branchless packed-tail mask was slower than the conditional form.
- Removing minimizer-position bookkeeping produced worse Phase 1 code on the
  benchmark CPU.
- Replacing the `k=31` packed k-mer initializer with an unaligned load and
  byte-swap was not a repeatable Phase 2 win (9.409 s candidate mean versus
  9.372 s checkpoint mean).
- Moving rare overflow histogram updates directly into the insertion branch
  made common-path code layout worse (9.798 s candidate mean versus 9.372 s
  checkpoint mean), despite removing thread-local flag accesses.
- Assigning the minimal count width to every individual KFF record reduced the
  one-pair file from 38.40 GB to 34.64 GB, but raised Phase 2 from 30.09 s to
  31.58–31.79 s. Batch-adaptive widths remain the faster policy.

The remaining difference is Phase 1 on the human pair: Tuna's 11.25 billion
superkmers exceed KMC's 8.41 billion. Future work should target boundary
density, packed parser-to-consumer transfer, and temporary-write volume without
reintroducing Phase 2 skew.

## Human3 completion gate

The final scale gate uses all 36 FASTQ files in
`/scratch4/rob/tuna_benchmarking/human3` (614,082,898,328 compressed bytes),
canonical `k=31`, 16 threads, and a 256 GB declared memory budget. Both tools
skip result serialization. The exact input manifest is
`benchmark/human3_all.fof`.

| Tool | Phase 1 | Phase 2 | Process wall | Peak RSS | Total | Distinct |
|---|---:|---:|---:|---:|---:|---:|
| Tuna optimization-round2 | 494.03 s | 306.50 s | 835.66 s | 8,753,464 KB | 516,924,379,564 | 20,868,636,896 |
| KMC 3 | 500.08 s | 572.19 s | 1,072.36 s | 249,614,536 KB | 516,924,379,564 | 20,868,636,896 |

Tuna is 22.1 percent faster overall, 1.2 percent faster in Phase 1, and 46.4
percent faster in Phase 2 while using 96.5 percent less peak memory. Total and
distinct counts match exactly.

The representative four-pair subset also exposes the full 16-thread topology:
Tuna disk mode completed in 201.05 s versus KMC's 259.53 s, with exact totals
of 122,567,140,288 and 6,985,961,483 distinct. A shared gzip
producer/consumer harness fixed the few-file in-memory scheduler; on one pair,
Tuna fell from 245.56 s to 87.87 s and beat KMC's 116.16 s. The four-pair
in-memory experiment reached 199,487,684 KB RSS and saved only another 12.3 s,
so the conservative disk selection remains appropriate for the full dataset.

Authoritative logs are retained under:

```text
/scratch4/rob/tuna_benchmarking/results/round2_human3_20260724/
```

## Human3 binary-output gate

The final apples-to-apples gate repeats the 614 GB, 36-file run with real
binary counts. Tuna writes standards-compliant adaptive-width KFF; KMC writes
its native binary database. Neither command uses Tuna `-co` or KMC `-w`.

| Tool/output | Phase 1 | Phase 2 + output | Process wall | Peak RSS | Output bytes |
|---|---:|---:|---:|---:|---:|
| Tuna adaptive KFF | 496.70 s | 294.80 s | 822.27 s | 13,222,728 KB | 212,223,183,420 |
| KMC native DB | 542.07 s | 604.55 s | 1,146.69 s | 250,576,008 KB | 204,998,649,924 |

Tuna is 28.3 percent faster overall, 8.4 percent faster in Phase 1, and 51.2
percent faster in Phase 2 including output, while using 94.7 percent less peak
memory. Its portable KFF is 3.5 percent larger than KMC's native database.
Both report exactly 516,924,379,564 total and 20,868,636,896 distinct k-mers.

An independent full KFF reader pass decoded every record and recovered those
same totals; the maximum encoded count was 54,689,425. On the one-pair subset,
Tuna adaptive KFF completed in 91.72 s versus 118.17 s for KMC native output
and 119.62 s for KMC KFF output. Tuna's 38,395,090,244-byte pair KFF was also
smaller than both KMC's 38,552,359,760-byte native database and its
46,181,052,101-byte KFF.

Authoritative real-output logs are retained under:

```text
/scratch4/rob/tuna_benchmarking/results/round2_output_20260724/
```

## Thread scaling and final continuation

The output-enabled `m=15` checkpoint was also measured at 8 and 32 threads
with the same 36-file manifest and 256 GB RAM budget:

| Threads | Tuna Phase 1 | Tuna Phase 2 + KFF | Tuna wall | KMC Phase 1 | KMC Phase 2 + DB | KMC wall | Tuna wall advantage |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 955.91 s | 595.02 s | 1,586.36 s | 992.32 s | 983.29 s | 1,975.66 s | 19.7% |
| 16 | 496.70 s | 294.80 s | 822.27 s | 542.07 s | 604.55 s | 1,146.69 s | 28.3% |
| 32 | 271.03 s | 260.14 s | 564.63 s | 463.38 s | 398.38 s | 861.82 s | 34.5% |

Tuna Phase 1 scales strongly through 32 threads. Phase 2 scales nearly linearly
from 8 to 16 threads but only improves from 294.80 to 260.14 seconds from 16
to 32, identifying counting and sustained output traffic as the remaining
scaling limit.

An output-enabled two-file minimizer sweep selected `m=17`, not the
count-only winner:

| `m` | Mode | Phase 1 | Phase 2 | Wall | Superkmers |
|---:|---|---:|---:|---:|---:|
| 13 | count only | 63.06 s | 40.93 s | 104.10 s | 3,570,084,037 |
| 15 | count only | 64.31 s | 26.03 s | 90.42 s | 3,914,225,849 |
| 17 | count only | 64.69 s | 17.72 s | 82.53 s | 4,349,968,508 |
| 19 | count only | 63.59 s | 17.24 s | 80.93 s | 4,906,939,931 |
| 21 | count only | 63.49 s | 19.00 s | 82.60 s | 5,644,207,915 |
| 17 | KFF output | 60.70 s | 21.54 s | 82.50 s | 4,349,968,508 |
| 19 | KFF output | 61.26 s | 23.53 s | 85.03 s | 4,906,939,931 |

All sweep rows report exactly 31,697,312,651 total and 3,848,420,222 distinct
k-mers. The shorter `m=13` reduces record count but loses partition entropy;
`m=19` narrowly wins count-only while its additional records make real output
slower. The production default is therefore `m=17`.

The retained rolling-hash, prefetch, KFF, cleanup, and `m=17` changes were then
run through the full 32-thread output gate:

| Tool/output | Phase 1 | Phase 2 + output | Process wall | Output bytes |
|---|---:|---:|---:|---:|
| Tuna final candidate KFF | 272.43 s | 198.25 s | 471.15 s | 212,496,712,212 |
| Tuna `m=15` checkpoint KFF | 271.03 s | 260.14 s | 564.63 s | 212,230,376,419 |
| KMC native DB | 463.38 s | 398.38 s | 861.82 s | 204,998,649,924 |

The final candidate is another 16.6 percent faster than the preceding Tuna
checkpoint and 45.3 percent faster than KMC overall. It reduces Tuna Phase 2
plus output by 23.8 percent without materially changing Phase 1. Tuna and KMC
again report exactly 516,924,379,564 total and 20,868,636,896 distinct
canonical k-mers. An independent KFF API scan decoded all 20,868,636,896
records, recovered the same 516,924,379,564 count sum, and observed a maximum
count of 54,689,425.

Continuation logs are retained under:

```text
/scratch4/rob/tuna_benchmarking/results/round2_output_scaling_20260724/
/scratch4/rob/tuna_benchmarking/results/round2_m_count_sweep2_20260724/
/scratch4/rob/tuna_benchmarking/results/round2_m_output_sweep_20260724/
/scratch4/rob/tuna_benchmarking/results/round2_final_candidate_20260724/
```
