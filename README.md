# tuna

> **Branch `feat/prot-support` — work in progress.**
> This branch adds `tuna-prot`, a protein k-mer counter built on the same minimizer-superkmer pipeline as tuna.
> It is not yet merged or released. Expect API and CLI changes.

**tuna** is a fast, streaming k-mer counter for FASTA/FASTQ input.
It partitions k-mers by minimizer into superkmer files, then counts them using a streaming hash table — keeping memory usage low and throughput high.

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

---

## How it works

tuna runs a two-phase pipeline:

1. **Partition (Phase 1)** — streams each input file through a minimizer iterator. Whenever the minimizer changes, the current superkmer is flushed to a per-partition binary file (on disk if unsufficient RAM budget). This groups k-mers that share a minimizer into the same bucket. The number of partitions is auto-tuned from input size (targeting ~2 MB input per partition) or set explicitly with `-n`.

2. **Count (Phase 2)** — replays each partition, upserting every k-mer into a Kache-hash table with increment semantics. Each partition is processed independently, so the hash table only ever holds one partition's k-mers at a time.

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

The default build embeds a dispatch table for **odd k in `[11, 31]`** — you can pass any of those values at runtime with no extra flags.

To use **k > 31**, **even k**, or any k up to 256, compile with explicit `-DFIXED_K` and `-DFIXED_M`:

```bash
cmake .. -DFIXED_K=63  -DFIXED_M=21
cmake .. -DFIXED_K=127 -DFIXED_M=21
cmake .. -DFIXED_K=256 -DFIXED_M=31
```

This produces a single-instantiation binary locked to that (k, m) pair. Those values become the defaults — `-k` and `-m` at runtime must match or the binary exits with an error. Build is faster and the binary is smaller. The same mechanism works for k ≤ 31 if you want a leaner binary:

```bash
cmake .. -DFIXED_K=31 -DFIXED_M=21
```

If you want to experiment with multiple (k, m) combinations, we recommend to use a separate build directory for each:

```bash
cmake -S . -B build_k31_m21  -DFIXED_K=31  -DFIXED_M=21  && cmake --build build_k31_m21  --target tuna -j$(nproc)
cmake -S . -B build_k63_m21  -DFIXED_K=63  -DFIXED_M=21  && cmake --build build_k63_m21  --target tuna -j$(nproc)
cmake -S . -B build_k127_m25 -DFIXED_K=127 -DFIXED_M=25  && cmake --build build_k127_m25 --target tuna -j$(nproc)
```

Each directory contains its own `tuna` binary: `build_k31_m21/tuna`, `build_k63_m21/tuna`, etc.

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
| `-k` | `<int>` | `31` | k-mer length. Odd values in `[11, 31]` in the default build; any value in `[2, 256]` when compiled with `-DFIXED_K=k` (must match compile-time value if set, default is compile-time value if set) |
| `-m` | `<int>` | `21` | Minimizer length. Any value in `[1, min(k-1, 32)]`. `m=21` is a good default; use `m=23`–`25` for highly repetitive or low-complexity data (e.g. individual human genomes). Must match `-DFIXED_M` if set, default is compile-time value if set |
| `-t` | `<int>` | `1` | Number of threads. Phase 1 parallelises over input files; Phase 2 over partitions |
| `-ci` | `<int>` | `1` | Minimum count to report |
| `-cx` | `<int>` | `max` | Maximum count to report |
| `-ram` | `<int>` | auto | RAM budget in GB. Controls whether the in-memory or disk pipeline is used, and sizes write buffers accordingly. Set lower than physical RAM to leave headroom for other processes, or higher to force the in-memory pipeline |
| `-w` | `<dir>` | next to output | Working directory for temporary partition files. |
| `-kff` | — | off | Write output in [KFF binary format](https://github.com/Kmer-File-Format/kff-reference) instead of TSV. Auto-detected from a `.kff` output extension. |
| `-h` / `--help` | — | — | Print usage |

<details>
<summary><strong>Advanced / benchmarking flags</strong></summary>

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `-n` | `<int>` | auto | Number of partitions. Auto-tuned to ~2 MB input/partition when omitted |
| `-hp` | — | off | Hide progress messages (phase timings are always emitted to stderr) |
| `-kt` | — | off | Keep temporary partition files after the run |
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

> **Large genomes** — counting a human-scale genome (3 Gbp) produces ~500 million unique k-mers. In TSV this reaches ~20–30 GB; in KFF binary (~12 bytes/k-mer) it is ~6 GB.

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

[K-mer File Format](https://github.com/Kmer-File-Format/kff-reference) binary output. Each k-mer is stored as a 2-bit packed sequence (A=0, C=1, G=2, T=3, MSB-first) with a 4-byte big-endian count. The file is marked `canonical=true` and `unique=true`. Roughly 3–4× smaller than TSV for k=31.

KFF files can be read with [kff-cpp-api](https://github.com/Kmer-File-Format/kff-cpp-api) or any other KFF-compatible tool.

---

Only k-mers with counts in `[ci, cx]` are written. The canonical (lexicographically smaller of forward/reverse-complement) form of each k-mer is reported.

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

Comparison with [KMC 3.2.4](https://github.com/refresh-bio/KMC), k=31, m=21, 8 threads, on a cluster node.
Each row shows the **median wall time** over per-file runs (100 files for bacteria/metagenomes, 10 for human and Tara).

| dataset | type | tuna median | KMC median | speedup | tuna p1 | tuna p2 |
|---------|------|-------------|------------|---------|---------|---------|
| *E. coli* | genomes (plain FASTA) | 0.50 s | 1.24 s | **2.5×** | 0.14 s | 0.30 s |
| *Salmonella* | genomes (gz) | 0.47 s | 1.25 s | **2.7×** | 0.12 s | 0.29 s |
| Gut | metagenome assemblies (plain FASTA) | 0.23 s | 0.75 s | **3.3×** | 0.07 s | 0.13 s |
| Human | genomes (gz) | 134 s | 208 s | **1.5×** | 42 s | 88 s |
| Tara | metagenome reads (gz, 5.9 GB) | 93 s | 177 s | **1.9×** | 28 s | 62 s |

tuna is consistently faster than KMC across all dataset types.
Memory usage scales with unique k-mers per partition rather than total input size.

![Per-file benchmark: wall time distributions, phase breakdown, and speedup across 5 datasets](benchmark/datasets.png)

---

## tuna-prot — protein k-mer counter (work in progress)

`tuna-prot` is a companion binary for counting k-mers in **protein sequences** (FASTA, plain or gzipped).
It shares the same minimizer-superkmer pipeline as tuna but uses a 5-bit amino acid encoding (k ≤ 12) and operates on the 20 standard amino acids.
Sequences are split on ambiguous residues (B, J, O, U, X, Z, `*`, `-`); multi-line FASTA is fully supported.

### Build

```bash
cd build && cmake --build . --target tuna_prot -j$(nproc)
```

The `tuna_prot` binary will be at `build/tuna_prot`.

### Usage

```
tuna_prot [options] <file1.fa> [file2.fa ...]
```

Input files may be plain FASTA or gzip FASTA (`.gz`).
Prefix a file with `@` to treat it as a file-of-filenames.

| Flag | Default | Description |
|------|---------|-------------|
| `-k <int>` | `7` | k-mer length `[4, 12]` |
| `-m <int>` | `5` | Minimizer length `[2, k-1]` |
| `-t <int>` | `4` | Number of threads |
| `-n <int>` | auto | Number of partitions (power of 2) |
| `-o <file>` | stdout | Output file (`-` = stdout) |
| `-w <dir>` | `.` | Work directory for temp files |
| `-ram <GB>` | auto | RAM budget |
| `-hp` | off | Hide progress |
| `-kt` | off | Keep temp partition files |
| `-tp` | off | Stop after phase 1 |

Output is TSV: `<kmer>\t<count>`, one k-mer per line.

### Example

```bash
tuna_prot -k 7 -m 4 -t 4 proteins.fasta > counts.tsv
```
