# tuna

**tuna** is a fast, streaming k-mer counter for FASTA/FASTQ input.
It partitions k-mers by minimizer into superkmer files, then counts them using a streaming hash table — keeping memory usage low and throughput high.

It uses [kache-hash](https://github.com/vicLeva/kache-hash) as its streaming k-mer hash table.
Phase 1 parsing uses a C++ port of [helicase](https://github.com/imartayan/helicase) (SIMD FASTA/FASTQ, ~5 GB/s), and minimizer hashing uses a C++ port of [simd-minimizers](https://github.com/Daniel-Liu-c0deb0t/simd-minimizers) (canonical ntHash, two-stack sliding window minimum).

---

## Table of contents

- [How it works](#how-it-works)
- [Partition strategies](#partition-strategies)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Output format](#output-format)
- [Benchmarks](#benchmarks)

---

## How it works

tuna runs a two-phase pipeline:

1. **Partition (Phase 1)** — streams each input file through a minimizer iterator. Whenever the minimizer changes, the current superkmer is flushed to a per-partition binary file on disk. This groups k-mers that share a minimizer into the same bucket. An optional **Phase 0** pre-scan can be run first to build a load-balanced partition table (see [Partition strategies](#partition-strategies)).

2. **Count (Phase 2)** — replays each partition through a `Kmer_Window`, upserting every k-mer into a `Streaming_Kmer_Hash_Table` with increment semantics. Each partition is processed independently, so the hash table only ever holds one partition's k-mers at a time.

3. **Output (Phase 2, cont.)** — iterates the table, applies `-ci`/`-cx` count filters, and writes `<kmer>\t<count>` to the output file.

Partitions are processed in parallel across threads (up to `-n` partitions at a time), keeping peak memory proportional to a single partition's k-mer set.

---

## Partition strategies

tuna supports three strategies for assigning k-mers to partitions, selected at runtime:

| Flag | Strategy | Pre-scan | Description |
|------|----------|----------|-------------|
| *(none / `-hash`)* | **Hash** (default) | No | `minimizer_hash % num_partitions`. Uses canonical ntHash (fwd XOR rev). No pre-scan; produces well-balanced partitions in practice. |
| `-kmtricks` | **kmtricks repart** | Yes | Partition scheme from [kmtricks](https://github.com/tlemane/kmtricks). Pre-scans input to count minimizer frequencies, then builds a load-balanced minimizer→partition table using LPT greedy assignment. |
| `-kmc` | **KMC signature** | Yes | Partition scheme from [KMC](https://github.com/refresh-bio/KMC). Uses KMC's norm-filtered canonical m-mer signature (lexicographic order + validity filter). Pre-scans input to build a load-balanced signature→partition table. |

All three strategies produce identical k-mer sets and counts.

---

## Dependencies

- C++20 compiler: [GCC](https://gcc.gnu.org/) >= 9.1 or [Clang](https://clang.llvm.org) >= 9.0
- [CMake](https://cmake.org/) >= 3.17
- [zlib](https://zlib.net/)

**Debian/Ubuntu:**
```bash
sudo apt-get install build-essential cmake zlib1g-dev
```

**Fedora/RHEL:**
```bash
sudo dnf install gcc-c++ cmake zlib-devel
```

**macOS:**
```bash
brew install llvm cmake zlib
```

---

## Installation

```bash
git clone https://github.com/vicLeva/tuna.git
cd tuna/
mkdir build && cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target tuna -j$(nproc)
```

> **Note:** Always pass `-DCMAKE_BUILD_TYPE=Release` explicitly. If a stale `CMakeCache.txt` exists in the build directory from a previous Debug configuration, plain `cmake ..` will reuse the cached build type and produce an unoptimized binary (~6× slower).

The `tuna` binary will be at `build/tuna`.

### Optional compile-time flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DINSTANCE_COUNT=N` | 32 | Number of k values instantiated at compile time |
| `-DFIXED_K=N` | — | Compile for a single k value only (smaller binary) |

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
| `-k` | `<int>` | `31` | k-mer length. Any odd value in `[11,31]` (fits in 64-bit word) |
| `-m` | `<int>` | `17` | Minimizer length. Any odd value in `[9, k-2]` (must be odd and < k) |
| `-n` | `<int>` | `32` | Number of partitions (superkmer buckets written to disk) |
| `-t` | `<int>` | `1` | Number of threads. Phase 1 parallelises over input files (capped at number of files); Phase 2 parallelises over partitions (capped at `-n`) |
| `-ci` | `<int>` | `1` | Minimum count to report |
| `-cx` | `<int>` | `max` | Maximum count to report |
| `-w` | `<dir>` | next to output | Working directory for temporary partition files |
| `-hash` | — | (default) | Hash partition strategy: `minimizer_hash % num_partitions`, no pre-scan |
| `-kmtricks` | — | off | kmtricks load-balanced partition strategy (requires Phase 0 pre-scan) |
| `-kmc` | — | off | KMC signature partition strategy (requires Phase 0 pre-scan) |
| `-s` | `<int>` | `9` | KMC signature length (`-kmc` only). Range: `[5, 11]` and must be < k |
| `-hp` | — | off | Hide progress messages (phase timings are always emitted to stderr) |
| `-kt` | — | off | Keep temporary partition files after the run (useful for benchmarking) |
| `-tp` | — | off | Stop after partitioning — Phase 1 only (useful for benchmarking partition speed) |
| `-h` / `--help` | — | — | Print usage |

### Examples

Count k-mers in a reference genome, k=31, 4 threads, report only k-mers seen at least twice:

```bash
tuna -k 31 -t 4 -ci 2 genome.fa counts.tsv
```

Count from a list of FASTA files using the kmtricks balanced partitioning:

```bash
tuna -k 31 -t 8 -kmtricks @genomes.list counts.tsv
```

Count from a list of FASTA files using the simple hash strategy (no pre-scan):

```bash
tuna -k 31 -t 8 -hash @genomes.list counts.tsv
```

---

## Output format

Plain text, tab-separated, one k-mer per line:

```
ACGTACGTACGTACGTACGTACGTACGTACG	42
TGCATGCATGCATGCATGCATGCATGCATGC	7
...
```

Only k-mers with counts in `[ci, cx]` are written. The canonical (lexicographically smaller of forward/reverse-complement) form of each k-mer is reported.

---

## Benchmarks

Comparison with [KMC 3.2.4](https://github.com/refresh-bio/KMC), k=31, 8 threads, `-kmc` partition strategy.

| dataset | tuna wall | tuna RSS | KMC wall | KMC RSS | time ratio | mem ratio |
|---------|-----------|----------|----------|---------|------------|-----------|
| *E. coli* ×10 (~47 MB) | 0:02 | 0.4 GB | 0:02 | 0.8 GB | 1.0× | **0.5×** |
| *E. coli* ×3682 (~17 GB) | 3:48 | 59.4 GB | 2:01 | 11.2 GB | 1.9× | 5.3× |
| human ×1 (~3 GB) | 2:12 | 10.1 GB | 2:20 | 11.2 GB | **0.9×** | **0.9×** |
| human ×10 (~30 GB) | 6:52 | 70.5 GB | 5:45 | 11.2 GB | 1.2× | 6.3× |

tuna's memory usage scales with input size (streaming approach — hash tables grow with the number of distinct k-mers per partition). KMC uses a roughly constant ~11 GB regardless of input size thanks to its on-disk intermediate representation. tuna is competitive on small-to-medium datasets and performs better than KMC on a single human genome.

Both tools produce identical k-mer sets and counts (verified on all three partition strategies).
