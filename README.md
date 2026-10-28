# tuna

**tuna** is a fast, streaming k-mer counter for FASTA/FASTQ input.
It partitions minimizer superkmers, then counts them with a partition-local
canonical rolling-hash table, keeping memory usage low and throughput high.

It uses [unordered_dense](https://github.com/martinus/unordered_dense) as its streaming k-mer hash table.
Phase 1 parsing uses a C++ port of [helicase](https://github.com/imartayan/helicase) (FASTX parser), and minimizers are computed with a two-stack sliding-window minimum over canonical m-mers, under a regular random minimizer scheme.

---

## Table of contents

- [How it works](#how-it-works)
- [Dependencies](#dependencies)
- [Installation](#installation)
  - [Choosing m](#choosing-m)
  - [k range and compile-time specialisation](#k-range-and-compile-time-specialisation)
- [Usage](#usage)
- [Output format](#output-format)
  - [TSV (default)](#tsv-default)
  - [KFF binary](#kff-binary--kff-or-kff-extension)
- [C++ library API](#c-library-api)
- [Benchmarks](#benchmarks)
  - [Datasets](#datasets)
  - [Reproducing](#reproducing)
  - [Results](#results)

---

## How it works

tuna runs a two-phase pipeline:

1. **Partition (Phase 1)** — overlaps gzip decompression, FASTX parsing, and minimizer iteration. Whenever the minimizer changes, the current packed superkmer is flushed to a binary partition in memory or on disk. tuna estimates the superkmer size up front and keeps phase 1 in memory when it fits the RAM budget; if the budget is reached while running, the partitions held in memory are spilled to disk and counting continues from there, so a run never exceeds its budget. Running with `-dbg` reports when a spill happens. The number of partitions is auto-tuned from input size or set explicitly with `-n`.

2. **Count (Phase 2)** — replays each partition and inserts its k-mers into a hash table, incrementing the count of each. Identical superkmers are collapsed first, so a repeated one is counted once and its multiplicity added in a single step. Each worker holds one partition table at a time.

3. **Output (Phase 2, cont.)** — iterates the table, applies `-ci`/`-cx` count filters, and writes results to the output file in TSV or [KFF](#output-format) format.

Partitions are processed in parallel across threads (up to `-n` partitions at a time), keeping peak memory proportional to a single partition's k-mer set.

---

## Dependencies

- **Platform: Linux or macOS** (x86\_64 and AArch64/Apple Silicon)
- C++20 compiler: [GCC](https://gcc.gnu.org/) >= 9.1 or [Clang](https://clang.llvm.org) >= 9.0
- [CMake](https://cmake.org/) >= 3.17

All other dependencies ([zlib-ng](https://github.com/zlib-ng/zlib-ng), [kff-cpp-api](https://github.com/Kmer-File-Format/kff-cpp-api)) are fetched and built automatically by CMake — no manual installation needed.

**Debian/Ubuntu:**
```bash
sudo apt-get install build-essential cmake
```

**Fedora/RHEL:**
```bash
sudo dnf install gcc-c++ cmake
```

**macOS:**
```bash
brew install llvm cmake
```

---

## Installation

```bash
git clone https://github.com/vicLeva/tuna.git
cd tuna/
mkdir build && cd build/
cmake ..
make -j$(nproc)
```

The `tuna` binary will be at `build/tuna`.

### Choosing m

You do not. Left out, `-m` is derived from `k`: **13 or less below `k=31`, 17 at or
above**. 

A build therefore instantiates one m per k. If for some reason (we do not recommend it) you would want to choose `-m` freely, add the following to the `cmake` command:

```bash
cmake .. -DTUNA_ALL_M=ON      # every odd m for every k
```

Expect the following `make` to be slow and memory-hungry: it will be 270 instantiations against 28.

### k range and compile-time specialisation

The default build accepts **odd k in `[11, 63]`, plus `k=127`**, each carrying
one minimizer length, `auto_m(k)`. That is 28 template instantiations.

**Most users want exactly this and should stop here.** Leave `-m` alone: the
value tuna picks is within a few percent of the best one, and overriding it
means building every m, which is ten times the instantiations and by far the
slowest build.

If the build is still too slow and you know your k in advance, pin it:

```bash
cmake -S . -B build_k31 -DFIXED_K=31 && cmake --build build_k31 --target tuna -j$(nproc)
```

One instantiation instead of 28, a smaller binary, and the fastest build.
The binary then accepts only `k=31`. Use a separate directory per k, which is
also how you keep several around at once: `build_k31/tuna`, `build_k63/tuna`.

#### Flags summary

| cmake flags | k accepted | m accepted | pairs |
|---|---|---|---|
| *(none)* | odd 11–63, 127 | `auto_m(k)` | 28 |
| `-DFIXED_K=K` | `K` | `auto_m(K)` | 1 |
| `-DTUNA_ALL_M=ON` | odd 11–63, 127 | any odd m | 270 |
| `-DFIXED_K=K -DTUNA_ALL_M=ON` | `K` | any odd m | ~12 |
| `-DFIXED_K=K -DFIXED_M=M` | **any k ≤ 256** | `M` | 1 |

Three notes on that table. The two switches do one thing each and are
independent: `FIXED_K` narrows which k, `TUNA_ALL_M` widens which m, so
sweeping m at a single k wants both.

The last row is the exception. The first four pick from a hardcoded list of k,
and `FIXED_K` can only remove entries from it, never add one — so
`-DFIXED_K=256` alone builds an empty table and fails at runtime. Giving
`FIXED_M` too switches to a branch with no table at all, which instantiates
just the pair you named and so reaches any `k <= 256`, even k included:

```bash
cmake .. -DFIXED_K=256 -DFIXED_M=31
```

`FIXED_M` does nothing without `FIXED_K`, and when both are set they win over
`TUNA_ALL_M`. Whatever the build, `-k` and `-m` at runtime must match what it
carries, or tuna exits saying which value it has.

<details>
<summary><strong>Other compile-time options</strong></summary>

**Debug build** — disables optimisations, enables debug symbols for gdb/valgrind:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

</details>

---

## Usage

```
tuna [options] <input1.fa [input2.fa ...]> <output_file>
tuna [options] @<input_list_file>          <output_file>
```

Input files can be FASTA or FASTQ, plain or gzipped.
Instead of listing files directly, you can pass `@list.txt` where `list.txt` is a newline-separated file of paths.

### Options

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `-k` | `<int>` | `31` | k-mer length. Odd values in `[11, 63]` plus `127` in the default build. Any value in `[2, 256]`, even ones included, needs a build with both `-DFIXED_K=k` and `-DFIXED_M=m`; `-k` must then match it. See [k range](#k-range-and-compile-time-specialisation). |
| `-t` | `<int>` | `1` | Number of threads. Compressed phase 1 pipelines decompression, parsing, and partitioning; phase 2 parallelizes over partitions. |
| `-ci` | `<int>` | `1` | Minimum count to report |
| `-cx` | `<int>` | `max` | Maximum count to report |
| `-ram` | `<int>` | auto | RAM budget in GB. Phase 1 starts in memory when the estimated superkmer size fits the budget and goes to disk otherwise; if the budget is reached mid-run, the partitions held in memory are spilled to disk and the run finishes on the disk path. Also sizes write buffers. Set lower than physical RAM to leave headroom for other processes, or higher to keep phase 1 in memory |
| `-w` | `<dir>` | next to output | Working directory for temporary partition files. |
| `-kff` | — | off | Write output in [KFF binary format](https://github.com/Kmer-File-Format/kff-reference) instead of TSV. Auto-detected from a `.kff` output extension. |
| `-b` | — | off | Disable canonical k-mers: count forward and reverse-complement strands independently. By default, a k-mer and its reverse complement are merged into a single count (the canonical, lexicographically smaller form is reported). Use `-b` when strand orientation matters or to match tools that count each strand separately. |
| `-h` / `--help` | — | — | Print usage |

<details>
<summary><strong>Advanced / benchmarking flags</strong></summary>

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `-m` | `<int>` | from `k` | Partition minimizer length, normally left alone: tuna derives it from `k` (**13 below `k=31`, 17 at or above**, clamped to `k-2`), within a few percent of the best value. A build carries only that one m per k, so another value needs `-DTUNA_ALL_M=ON` or `-DFIXED_M` — see [Choosing m](#choosing-m). Odd values only. |
| `-n` | `<int>` | auto | Number of partitions. Auto-tuned to ~2 MB input/partition when omitted |
| `-hp` | — | off | Hide progress messages (phase timings are always emitted to stderr) |
| `-kt` | — | off | Keep temporary partition files after the run |
| `-co` | — | off | Count only: skip k-mer serialization while still reporting total and distinct k-mer counts |
| `-tp` | — | off | Stop after partitioning — Phase 1 only |
| `-dbg` | — | off | Per-partition table summary + minimizer coverage CSV written to `<work_dir>/debug_min_coverage.csv` |

</details>

### Examples

Count k-mers in a reference genome, k=31, 4 threads:

```bash
tuna -k 31 -t 4 genome.fa counts.tsv
```

Count only k-mers seen at least twice:

```bash
tuna -k 31 -t 4 -ci 2 genome.fa counts.tsv
```

Count from a list of files:

```bash
tuna -k 31 -t 8 @genomes.list counts.tsv
```

Write KFF binary output (auto-detected from extension):

```bash
tuna -k 31 -t 8 @genomes.list counts.kff
```

Benchmark counting without serializing k-mers:

```bash
tuna -k 31 -t 8 -co @genomes.list /dev/null
```

> **Large genomes** — counting one human assembly (3.0 Gbp) yields about 2.5 billion distinct 31-mers. In TSV that is roughly 85 GB, so prefer KFF, or `/dev/null` with `-co` when benchmarking; at k=31, KFF uses 8 sequence bytes plus 1–4 count bytes per k-mer.

---

## Output format

### TSV (default)

Plain text, tab-separated, one k-mer per line:

```
ACGTACGTACGTACGTACGTACGTACGTACG	42
TGCATGCATGCATGCATGCATGCATGCATGC	7
...
```

### KFF binary (`-kff` or `.kff` extension)

[K-mer File Format](https://github.com/Kmer-File-Format/kff-reference) binary output. Each k-mer is stored as a 2-bit packed sequence (A=0, C=1, G=2, T=3) with its smallest lossless 1–4 byte big-endian count. Workers batch records by count width before writing KFF raw sections, avoiding per-record section changes. The file is marked `canonical=true` and `unique=true` (or `canonical=false` when `-b` is used). Roughly 3–4× smaller than TSV for k=31.

KFF files can be read with [kff-cpp-api](https://github.com/Kmer-File-Format/kff-cpp-api) or any other KFF-compatible tool.

---

Only k-mers with counts in `[ci, cx]` are written. By default, the canonical (lexicographically smaller of forward/reverse-complement) form of each k-mer is reported. With `-b`, the observed strand form is reported instead.

---

## C++ library API

tuna can be embedded directly in a C++ project

```cpp
#include <tuna/tuna.hpp>

// Collect all k-mers into a map (simple)
auto kmers = tuna::count_to<31>({"genome.fa"});   // std::unordered_map<std::string, uint32_t>

// Stream k-mers through a callback (memory-efficient)
tuna::count<31>({"genome.fa"}, [](std::string_view kmer, uint32_t count) {
    // called for every canonical k-mer; may run from multiple threads
});

// Large k: any value in [2, 256], both k and m are template parameters
tuna::count<127, 21>({"genome.fa"}, [](std::string_view kmer, uint32_t count) { ... });

// Options mirror the CLI flags; m defaults to 17, as on the command line
tuna::Options opts;
opts.threads   = 8;       // -t
opts.ram_gb    = 64;      // -ram, 0 = auto-detect
opts.canonical = false;   // -b: count each strand separately

tuna::count<31>({"reads.fq.gz"}, [](std::string_view kmer, uint32_t count) { ... }, opts);
```

The callback runs the same counting path as the CLI and differs only in where
counted k-mers go, so an API run and an equivalent `tuna` invocation produce
identical counts.

CMake integration:
```cmake
add_subdirectory(tuna)                                    # or use FetchContent
target_link_libraries(my_target PRIVATE tuna::tuna)
```

For a full walkthrough: CMake setup, FetchContent, container customisation, thread safety, see the **[wiki: Using tuna as a library](https://github.com/vicLeva/tuna/wiki/Using%E2%80%90tuna%E2%80%90as%E2%80%90a%E2%80%90library)**.

---

## Benchmarks

All numbers below were measured at `k=31`, 8 threads, with a 256 GB memory
budget, on a GenOuest node (two 32-core Xeon Gold 6338N at 2.20 GHz, 48 MB of
L3 per socket, 1 TB of RAM, Debian 12), against
[KMC 3.2.4](https://github.com/refresh-bio/KMC) and
[FastK](https://github.com/thegenemyers/FASTK). Every tool was asked to report
all k-mers from a count of one upwards, so the three produce the same set of
counts, and every run writes the tool's own binary format. tuna chose `m`
itself, which is 17 at `k=31`.

### Datasets

| dataset | organism | type | files | where to get it |
|---|---|---|---|---|
| Ecoli | *E. coli* | assemblies (plain FASTA) | 3500 | [Zenodo 6577997](https://zenodo.org/records/6577997) |
| Salmonella | *S. enterica* | assemblies (gz) | 25000 | [ENA2018-bacteria-661k](http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/) |
| Gut | gut MAGs | assemblies (gz) | 10000 | [HumGut](https://arken.nmbu.no/~larssn/humgut/) |
| Human | *H. sapiens* | assemblies (gz) | 60 | [HPP Year 1 Assemblies](https://github.com/human-pangenomics/HPP_Year1_Assemblies) |
| Tara | Tara Oceans (sea water) | reads (gz) | 12 | 6 ENA runs from [PRJEB4352](https://www.ebi.ac.uk/ena/browser/view/PRJEB4352), listed below |
| Gallus | *G. gallus* | reads (gz) | 12 | ENA `SRR105788`, `SRR105789`, `SRR105792`, `SRR105794`, `SRR197985`, `SRR197986` (paired) |
| Human3 | *H. sapiens* | reads (gz) | 36 | ENA `ERR174324`-`ERR174341` (paired) |
| HumanR | *H. sapiens* | reads (gz) | 1 | ENA `SRR622461`, forward reads only |

The *files* column is how many files each experiment counts, which for Ecoli,
Salmonella and Gut is a subset of a much larger public collection. The per-file
experiment uses the first 100 files of Ecoli, Salmonella and Gut and the first
10 of Human and Tara, and reports a median over them. The per-collection
experiment counts every file of the column in a single run. Gallus and Human3
are only ever counted whole. HumanR drives the coverage sweep.

The Tara files are surface-water samples, size fraction QQSS, from stations 7,
23, 30 and 11 of the Tara Oceans expedition. They come from six sequencing runs
of four samples. The ten files named below are the ones the per-file
experiment uses, the collection run taking twelve from the same six runs:

| ENA run | sample | station | files used |
|---|---|---|---|
| [`ERR315827`](https://www.ebi.ac.uk/ena/browser/view/ERR315827) | ERS327852 | 7 SUR1 | `AHX_AAGOSU_6_1_814P9ABXX`, `_6_2_814P9ABXX` |
| [`ERR318594`](https://www.ebi.ac.uk/ena/browser/view/ERR318594) | ERS329507 | 23 SUR2 | `AHX_AASOSU_1_2_62FGDAAXX` |
| [`ERR318583`](https://www.ebi.ac.uk/ena/browser/view/ERR318583) | ERS329507 | 23 SUR2 | `AHX_AASOSU_2_1_62FGDAAXX`, `_2_2_62FGDAAXX` |
| [`ERR318616`](https://www.ebi.ac.uk/ena/browser/view/ERR318616) | ERS329505 | 30 SUR2 | `AHX_ABEOSU_1_1_62J5HAAXX`, `_1_2_62J5HAAXX` |
| [`ERR318585`](https://www.ebi.ac.uk/ena/browser/view/ERR318585) | ERS329505 | 30 SUR2 | `AHX_ABEOSU_2_1_62J5HAAXX`, `_2_2_62J5HAAXX` |
| [`ERR1726642`](https://www.ebi.ac.uk/ena/browser/view/ERR1726642) | ERS488262 | 11 | `AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean` |

The submitted files can be fetched directly:

```
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR315/ERR315827/AHX_AAGOSU_6_1_814P9ABXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR315/ERR315827/AHX_AAGOSU_6_2_814P9ABXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318594/AHX_AASOSU_1_2_62FGDAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318583/AHX_AASOSU_2_1_62FGDAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318583/AHX_AASOSU_2_2_62FGDAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318616/AHX_ABEOSU_1_1_62J5HAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318616/AHX_ABEOSU_1_2_62J5HAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318585/AHX_ABEOSU_2_1_62J5HAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR318/ERR318585/AHX_ABEOSU_2_2_62J5HAAXX.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR172/ERR1726642/AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq.gz
```

### Reproducing

Every script that produced these numbers is in [`scripts/`](scripts/), one per
experiment, documented in [`scripts/README.md`](scripts/README.md). They share
`bench_common.sh` and are resumable: anything already in the output CSV is
skipped, so an interrupted run continues where it stopped.

### Results

Median wall time per file, each file counted in its own run, each tool writing
its own compact format (KFF for tuna, native databases for the others).

| dataset | tuna | KMC3 | FastK |
|---|---|---|---|
| Ecoli | **0.22 s** | 0.58 s | 0.37 s |
| Salmonella | **0.21 s** | 0.59 s | 0.41 s |
| Gut | **0.12 s** | 0.44 s | 0.23 s |
| Human | **27.5 s** | 33.8 s | 104.5 s |
| Tara | **33.4 s** | 37.6 s | 127.3 s |

The same collections counted in a single run over all of their files, which is
what counting a collection rather than its samples costs:

| dataset | files | tuna | KMC3 | FastK |
|---|---|---|---|---|
| Ecoli | 3500 | **29.6 s** | 127.4 s | 62.9 s |
| Salmonella | 25000 | **5.3 min** | 15.9 min | failed |
| Gut | 10000 | **1.8 min** | 2.8 min | 4.8 min |
| Human | 60 | **7.2 min** | 24.6 min | 11.2 min |
| Tara | 12 | 8.0 min | **5.9 min** | timeout |

KMC fills the budget it is given on every one of them, about 238 GB, where tuna
peaks between 4.8 and 62 GB and writes nothing at all to disk on Ecoli, Gut and
Tara, counting them from end to end in memory.

Tara is the one collection where sorting wins. Its 12 read sets hold 22.5
billion distinct k-mers at an abundance of 1.6, so the tables stay out of cache
and there is next to nothing for deduplication to collapse, which is exactly
the case where a sort-based counter is expected to come out ahead.

The two largest collections, counted whole:

| dataset | input | tuna | KMC3 | FastK | tuna RSS | KMC3 RSS | FastK RSS | tuna disk | KMC3 disk | FastK disk |
|---|---|---|---|---|---|---|---|---|---|---|
| Gallus | 24.0 GB | **2.5 min** | 3.5 min | timeout | **34.9 GB** | 238.5 GB | - | **0 GB** | 24.2 GB | - |
| Human3 | 614 GB | **59.0 min** | 67.1 min | 71.2 min | **3.8 GB** | 238.9 GB | 192.1 GB | 707 GB | 523 GB | 510 GB |

Both FastK runs on Gallus exceeded the six-hour limit. Gallus fits the memory
budget, so tuna counts it without touching the disk at all. Human3 does not:
tuna switches to its disk pipeline and counts 20.9 billion distinct k-mers in
3.8 GB of RAM, trading that for the temporary data in the last columns, which
all three tools write by the hundreds of gigabytes.
