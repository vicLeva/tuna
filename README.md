# tuna

**tuna** is a fast, streaming k-mer counter for FASTA/FASTQ input.
It partitions minimizer superkmers, then counts them with a partition-local
canonical rolling-hash table, keeping memory usage low and throughput high.

It uses [kache-hash](https://github.com/jamshed/kache-hash) as its streaming k-mer hash table.
Phase 1 parsing uses a C++ port of [helicase](https://github.com/imartayan/helicase) (SIMD FASTX parser), and minimizer hashing uses a C++ port of [simd-minimizers](https://github.com/rust-seq/simd-minimizers) (canonical ntHash, two-stack sliding window minimum).

---

## Table of contents

- [How it works](#how-it-works)
- [Dependencies](#dependencies)
- [Installation](#installation)
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

1. **Partition (Phase 1)** — overlaps gzip decompression, FASTX parsing, and minimizer iteration. Whenever the minimizer changes, the current packed superkmer is flushed to a binary partition in memory or on disk. The number of partitions is auto-tuned from input size or set explicitly with `-n`.

2. **Count (Phase 2)** — replays each partition, routing canonical rolling k-mer hashes into a Kache-hash table with increment semantics. This is independent of the shorter partition minimizer, avoiding minimizer-induced bucket skew. Each worker holds one partition table at a time.

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

### k range and compile-time specialisation

The default build embeds a dispatch table for **odd k in `[11, 63]`, plus
`k=127`**. You can pass any of those values at runtime with no extra flags.

To use **even k**, an odd k not in that table, or any k up to 256, compile with
explicit `-DFIXED_K` and `-DFIXED_M`:

```bash
cmake .. -DFIXED_K=63  -DFIXED_M=21
cmake .. -DFIXED_K=127 -DFIXED_M=21
cmake .. -DFIXED_K=256 -DFIXED_M=31
```

This produces a single-instantiation binary locked to that (k, m) pair. Those values become the defaults — `-k` and `-m` at runtime must match or the binary exits with an error. Build is faster and the binary is smaller. The same mechanism works for k ≤ 31 if you want a leaner binary:

```bash
cmake .. -DFIXED_K=31 -DFIXED_M=17
```

If you want to experiment with multiple (k, m) combinations, we recommend to use a separate build directory for each:

```bash
cmake -S . -B build_k31_m17  -DFIXED_K=31  -DFIXED_M=17  && cmake --build build_k31_m17  --target tuna -j$(nproc)
cmake -S . -B build_k63_m21  -DFIXED_K=63  -DFIXED_M=21  && cmake --build build_k63_m21  --target tuna -j$(nproc)
cmake -S . -B build_k127_m25 -DFIXED_K=127 -DFIXED_M=25  && cmake --build build_k127_m25 --target tuna -j$(nproc)
```

Each directory contains its own `tuna` binary: `build_k31_m17/tuna`, `build_k63_m21/tuna`, etc.

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
| `-k` | `<int>` | `31` | k-mer length. Odd values in `[11, 63]` plus `127` in the default build; any value in `[2, 256]` when compiled with `-DFIXED_K=k` (must match the compile-time value when set). |
| `-m` | `<int>` | `17` | Partition minimizer length. Any value in `[1, min(k-1, 32)]`; `17` balances compact superkmers with partition entropy for the default `k=31`. Phase-2 hash routing is independent of this value. Must match `-DFIXED_M` in a specialized build. |
| `-t` | `<int>` | `1` | Number of threads. Compressed phase 1 pipelines decompression, parsing, and partitioning; phase 2 parallelizes over partitions. |
| `-ci` | `<int>` | `1` | Minimum count to report |
| `-cx` | `<int>` | `max` | Maximum count to report |
| `-ram` | `<int>` | auto | RAM budget in GB. Controls whether the in-memory or disk pipeline is used, and sizes write buffers accordingly. Set lower than physical RAM to leave headroom for other processes, or higher to force the in-memory pipeline |
| `-w` | `<dir>` | next to output | Working directory for temporary partition files. |
| `-kff` | — | off | Write output in [KFF binary format](https://github.com/Kmer-File-Format/kff-reference) instead of TSV. Auto-detected from a `.kff` output extension. |
| `-b` | — | off | Disable canonical k-mers: count forward and reverse-complement strands independently. By default, a k-mer and its reverse complement are merged into a single count (the canonical, lexicographically smaller form is reported). Use `-b` when strand orientation matters or to match tools that count each strand separately. |
| `-h` / `--help` | — | — | Print usage |

<details>
<summary><strong>Advanced / benchmarking flags</strong></summary>

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
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

> **Large genomes** — counting a human-scale genome (3 Gbp) produces ~500 million unique k-mers. In TSV this reaches ~20–30 GB; at k=31, KFF uses 8 sequence bytes plus 1–4 count bytes per k-mer.

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
```

CMake integration:
```cmake
add_subdirectory(tuna)                                    # or use FetchContent
target_link_libraries(my_target PRIVATE tuna::tuna)
```

For a full walkthrough: CMake setup, FetchContent, container customisation, thread safety, see the **[wiki: Using tuna as a library](https://github.com/vicLeva/tuna/wiki/Using%E2%80%90tuna%E2%80%90as%E2%80%90a%E2%80%90library)**.

---

## Benchmarks

All numbers below were measured at `k=31`, `m=21`, 8 threads, with a 256 GB
memory budget, on a GenOuest node (4x8 Xeon E5-2660 at 2.20 GHz, 1.5 TB RAM,
CentOS 7), against [KMC 3.2.4](https://github.com/refresh-bio/KMC) and
[FastK](https://github.com/thegenemyers/FASTK). Every tool was asked to report
all k-mers from a count of one upwards, so the three produce the same set of
counts.

### Datasets

| dataset | organism | type | files | where to get it |
|---|---|---|---|---|
| Ecoli | *E. coli* | assemblies (plain FASTA) | 3682 | [Zenodo 6577997](https://zenodo.org/records/6577997) |
| Salmonella | *S. enterica* | assemblies (gz) | 10000 | [ENA2018-bacteria-661k](http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/) |
| Gut | gut MAGs | assemblies (gz) | 10000 | [HumGut](https://arken.nmbu.no/~larssn/humgut/) |
| Human | *H. sapiens* | assemblies (gz) | 60 | [HPP Year 1 Assemblies](https://github.com/human-pangenomics/HPP_Year1_Assemblies) |
| Tara | Tara Oceans (sea water) | reads (gz) | 10 | 6 ENA runs from [PRJEB4352](https://www.ebi.ac.uk/ena/browser/view/PRJEB4352), listed below |
| Gallus | *G. gallus* | reads (gz) | 12 | ENA `SRR105788`, `SRR105789`, `SRR105792`, `SRR105794`, `SRR197985`, `SRR197986` (paired) |
| Human3 | *H. sapiens* | reads (gz) | 36 | ENA `ERR174324`-`ERR174341` (paired) |
| HumanR | *H. sapiens* | reads (gz) | 1 | ENA `SRR622461`, forward reads only |

The per-file experiments use the first 100 files of Ecoli, Salmonella and Gut,
and the first 10 of Human and Tara. Gallus and Human3 are counted whole, in a
single run over the entire collection. HumanR drives the coverage sweep.

The Tara files are surface-water samples, size fraction QQSS, from stations 7,
23, 30 and 11 of the Tara Oceans expedition. They come from six sequencing runs
of four samples, of which the benchmark uses ten files:

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

Median wall time per file, each file counted in its own run. *binary* is each
tool writing its own compact format (KFF for tuna, native databases for the
others), *ASCII* is the whole cost of obtaining a plain-text table, including
the separate conversion pass KMC and FastK each require.

| dataset | tuna (bin) | KMC3 (bin) | FastK (bin) | tuna (ASCII) | KMC3 (ASCII) | FastK (ASCII) |
|---|---|---|---|---|---|---|
| Ecoli | **0.29 s** | 0.57 s | 0.44 s | **0.41 s** | 1.22 s | 1.11 s |
| Salmonella | **0.23 s** | 0.59 s | 0.42 s | **0.34 s** | 1.22 s | 1.08 s |
| Gut | **0.14 s** | 0.44 s | 0.26 s | **0.20 s** | 0.74 s | 0.56 s |
| Human | **31.2 s** | 39.3 s | 105.6 s | **77.9 s** | 229.8 s | failed |
| Tara | **35.1 s** | 40.6 s | 128.6 s | **76.7 s** | 192.8 s | 289.0 s |

FastK's `Tabex` conversion step segfaults on every Human file, so it has no
ASCII result there.

Whole collections counted in a single run, binary output:

| dataset | input | tuna | KMC3 | FastK | tuna RSS | KMC3 RSS | FastK RSS |
|---|---|---|---|---|---|---|---|
| Gallus | 22.3 GB | **2.7 min** | 3.5 min | timeout | **43 GB** | 238 GB | - |
| Human3 | 572 GB | **60.0 min** | 66.4 min | 71.2 min | **3 GB** | 239 GB | 192 GB |

Both FastK runs on Gallus exceeded the six-hour limit. On Human3 the input no
longer fits the memory budget, so tuna switches to its disk pipeline and counts
20.9 billion distinct k-mers in 3 GB of RAM, while still finishing first.
