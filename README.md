# tuna

Streaming k-mer counter with multiple partitioning strategies.

tuna reads FASTA/FASTQ files (plain or gzip-compressed), partitions k-mers
into temporary superkmer files, then counts them in parallel using
[kache-hash](https://github.com/COMBINE-lab/kache-hash).
Output is a TSV of `kmer<TAB>count` pairs.

## Features

- **Three partitioning strategies**: KMC-style signature (default), kmtricks-style minimizer repart, or simple hash
- **Parallel**: configurable thread count for both partitioning and counting phases
- **Low memory**: streaming design; RAM scales with hash table size, not input
- **Format support**: FASTA (`.fa`, `.fna`, `.fasta`) and FASTQ (`.fq`, `.fastq`), plain or `.gz`

## Build

Requires: CMake ≥ 3.17, C++20 compiler, zlib.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target tuna -j$(nproc)
```

The binary is at `build/tuna`.

> **Note**: always build in Release mode. A Debug build is ~6× slower.

### Optional compile-time flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DINSTANCE_COUNT=N` | 32 | Number of k values instantiated at compile time |
| `-DFIXED_K=N` | — | Compile for a single k value only (smaller binary) |

## Usage

```
tuna [options] <input> <output.tsv>
```

`<input>` is either a single file or `@filelist` where `filelist` is a
text file with one path per line.

### Options

| Flag | Description |
|------|-------------|
| `-k INT` | k-mer length (default: 31) |
| `-m INT` | minimizer length (default: 17) |
| `-n INT` | number of partitions (default: 32) |
| `-t INT` | thread count (default: 4) |
| `-ci INT` | minimum count threshold (default: 1) |
| `-cx INT` | maximum count threshold (default: unlimited) |
| `-w DIR` | working directory for temp files (default: system tmp) |
| `-kmc` | use KMC signature partitioning (default) |
| `-s INT` | KMC signature length, 5–11 (default: 9) |
| `-kmtricks` | use kmtricks-style minimizer repart |
| `-hash` | use simple `hash % n_partitions` partitioning |
| `-hp` | hide progress output |
| `-kt` | keep temp partition files after run |
| `-tp` | stop after partitioning (phase 1 only) |

### Examples

```bash
# Count k-mers in a single file
tuna -k 31 -t 8 reads.fa.gz counts.tsv

# Count from a list of files
tuna -k 31 -t 8 @filelist.txt counts.tsv

# Use kmtricks partitioning with 16 partitions
tuna -k 31 -n 16 -t 8 -kmtricks @filelist.txt counts.tsv
```

## Partitioning strategies

| Flag | Strategy | Notes |
|------|----------|-------|
| (none) / `-kmc` | **KMC signature** (default) | Pre-scan → LPT-balanced table; best load balance |
| `-kmtricks` | **kmtricks repart** | Pre-scan → LPT greedy table (REPART_BITS=10) |
| `-hash` | **Simple hash** | `hash % n_partitions`; no pre-scan, fastest phase 1 |

## Output format

Tab-separated, one k-mer per line:

```
ACGTACGT...    3
TGCATGCA...    1
```

## Benchmarking

Scripts in `benchmark/`:

```bash
# Compare tuna vs KMC
bash benchmark/bench.sh @filelist.txt results/

# Compare all three partitioning strategies across thread counts
bash benchmark/bench_part_modes.sh @filelist.txt results/ \
    K=31 THREADS="1 2 4 8" PARTS=32 REPS=3

# Plot results
python3 benchmark/plot_part_modes.py results/
```

Path overrides (for cluster use):

```bash
export TUNA_BIN=/path/to/tuna
export KMC_BIN=/path/to/kmc
export GNU_TIME=/usr/bin/time
bash benchmark/bench_part_modes.sh ...
```

## Dependencies

- [kache-hash](https://github.com/COMBINE-lab/kache-hash) — concurrent hash table (included in `include/kache-hash/`)
- [xxHash](https://github.com/Cyan4973/xxHash) — hashing (included in `external/xxHash/`)
- [wyhash](https://github.com/wangyi-fudan/wyhash) — hashing (included in `include/wyhash/`)
- zlib — gzip decompression (system library)
