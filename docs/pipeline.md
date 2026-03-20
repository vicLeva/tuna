# tuna — Data Flow: From Genome File to k-mer Counts

Diagram: `docs/pipeline_diagram.png` / `pipeline_diagram.pdf`

---

## Overview

tuna operates in two sequential phases:

- **Phase 1 — Partition**: stream input sequences through a minimizer iterator, cut them into *superkmers*, and write binary superkmer records to per-partition files on disk.
- **Phase 2 — Count + Write**: replay each partition file through a sliding k-mer window, insert every k-mer into a private hash table with increment semantics, then filter and write the result as a TSV.

Both phases are fully parallel across threads with no shared mutable state between partitions.

---

## Step 0 — Startup (`src/main.cpp`, `src/CLI.cpp`)

- CLI flags are parsed into a `Config` struct: `k`, `l` (minimizer length), `-n` (number of partitions), `-t` (threads), `-ci`/`-cx` (min/max count filters), `-w` (work directory), etc.
- **Auto-tune partitions** (when `-n` is omitted):
  - Target ≈ 2 MB of uncompressed sequence per partition.
  - `n = next_pow2(total_input_bytes / 2 MB)`, clamped to `[16, fd_limit − 32]`.
  - `.gz` files are multiplied by 6 to account for compression.
  - For multi-file runs with small average file size (< 50 MB), `n` is further capped at `next_pow2(n_files)` to avoid over-partitioning redundant genomes.
- **Template dispatch**: `dispatch(k, l, cfg)` resolves at runtime to a concrete `run<K,L>(cfg)` instantiation — all inner loops are monomorphised at compile time with no virtual dispatch overhead.

---

## Step 1 — Sequence Reading (`src/seq_source.hpp`, `src/fast_fasta.hpp`)

Two paths depending on file format:

| Format | Reader | Mechanism |
|---|---|---|
| Plain FASTA / FASTQ | `helicase` mmap | OS maps the file read-only; helicase delivers ACTG-only chunks directly (zero-copy) |
| `.gz` FASTA / FASTQ | `SeqReader` + 64 MB sliding window | `gzread` decompresses into the window; `split_actg()` excises non-ACTG characters (N, newlines, headers) |

**Output of this step**: a stream of **ACTG-only chunks** — contiguous DNA strings with no ambiguity codes, no headers, and no whitespace.

**Threading model**:
- Multiple files → threads steal files atomically (file-level work-stealing). Each thread runs its own helicase/SeqReader instance with no shared state.
- Single `.gz` file → 1 producer thread decompresses + splits; `n_threads − 1` consumer threads process batches. A bounded queue (capacity 32 batches) provides backpressure.

---

## Step 2 — Minimizer Window (`include/minimizer_window.hpp`, `include/nt_hash.hpp`)

For each ACTG chunk, a `MinimizerWindow<k,l>` runs a streaming **canonical l-mer minimizer** over the k-mer window of width `w = k − l` l-mers.

### Rolling hash — `Roller<l>` (ntHash)

- Maintains two hashes: `fwd_` (forward strand) and `rev_` (reverse-complement strand).
- `canonical() = fwd_ XOR rev_` — strand-independent; same hash for a k-mer and its reverse complement.
- **Init**: O(l) — seed both hashes from the first l bases.
- **Roll**: O(1) — XOR out the leaving base, XOR in the arriving base, with a cyclic bit rotation.

### Two-stack sliding minimum

- A prefix array `M_pre[]` and a running suffix minimum `M_suf` together answer *"what is the minimum canonical hash in the current window of w l-mers?"* in **O(1) per step**.
- `hash()` — O(1): returns `min(M_pre[pivot], M_suf)`.
- `min_lmer_pos()` — O(1): returns the absolute position of the minimizer l-mer within the sequence, tracked as a side-effect of the sliding minimum. This replaces the old O(superkmer_len) position-rescan that consumed ~22% of Phase 1 cycles.

---

## Step 3 — Superkmer Extraction (`src/partition_hash.hpp` → `extract_superkmers_from_actg`)

A **superkmer** is a maximal run of consecutive k-mers that all share the same minimizer hash value. Crucially, this means every k-mer inside a superkmer maps to the **same partition** and shares the **same minimizer hash** — which only needs to be computed once in Phase 2.

```
Initialise window on the first k bases.
prev_hash = canonical ntHash of the minimizer l-mer
prev_min_pos = absolute position of that l-mer

For each new incoming base:
    Advance window → new_hash, new_min_pos

    If new_hash ≠ prev_hash  OR  current superkmer length ≥ 255:
        Flush superkmer [sk_start … sk_end) to partition (prev_hash % n_parts)
        min_pos = prev_min_pos − sk_start       ← relative offset within superkmer
        kmer_count += sk_len − k + 1

        sk_start = pos − (k − 1)                ← overlap k−1 bases with next superkmer
        prev_hash    = new_hash
        prev_min_pos = new_min_pos
```

- The **k−1 base overlap** at superkmer boundaries ensures that the k-mers spanning a minimizer boundary are not lost.
- Superkmers are capped at 255 bases (uint8_t header field); the window flush also triggers at this limit.
- Typical superkmer length: `2k − l` bases (the range over which a minimizer l-mer can be the window minimum).

Each flushed superkmer is serialised by `SuperkmerWriter` and buffered in a per-thread, per-partition buffer (flush threshold = `max(4 KB, 64 MB / n_parts)`), then written to the shared partition file under a per-partition mutex.

---

## Step 4 — Superkmer Binary Format (`src/superkmer_io.hpp`)

```
┌──────────┬──────────┬────────────────────────────────────┐
│  1 byte  │  1 byte  │  ⌈len_bases / 4⌉ bytes             │
│ len_bases│  min_pos │  packed DNA  (4 bases/byte, MSB)   │
└──────────┴──────────┴────────────────────────────────────┘
```

- `len_bases`: number of bases in the superkmer (max 255).
- `min_pos`: 0-indexed start of the minimizer l-mer within the superkmer. Stored so Phase 2 can recompute the minimizer hash in **O(l)** instead of O(k), without re-running the full MinimizerWindow.
- **Packing**: kache encoding (A=0, C=1, G=2, T=3), 2 bits per base, most-significant-bit first. 4× smaller than ASCII. Phase 2 unpacks directly to kache-hash's internal base type — no ASCII round-trip.

One file `work_dir/hash_P.superkmers` is produced per partition `P`.

---

## Step 5 — Phase 2 Setup (`src/count.hpp` → `count_and_write`)

Threads are assigned partitions by striping: thread `t` processes partitions `t, t + n_threads, t + 2·n_threads, …` Each thread owns its partitions **exclusively** — no inter-partition locking, no atomics on hash table operations.

**Hash table initial size** per partition:

| Condition | Formula | Rationale |
|---|---|---|
| `n_files ≤ 10` | `clamp(2 × total_kmers / n_parts, 256K, 4M)` | Few files → total k-mers ≈ unique k-mers; dynamic estimate avoids oversizing |
| `n_files > 10` | `clamp(128M / n_parts, 256K, 4M)` | Many redundant files → total ≫ unique; use conservative fixed formula |

---

## Step 6 — Partition File Reading (`src/superkmer_io.hpp` → `SuperkmerReader`)

Each partition file is opened with:
```
open() → fstat() → mmap(MAP_PRIVATE | MAP_POPULATE) → madvise(SEQUENTIAL | WILLNEED)
```

- `MAP_POPULATE` prefaults all pages before the counting loop, so the OS loads the data into RAM upfront.
- `MADV_SEQUENTIAL | MADV_WILLNEED` hints to the kernel to prefetch ahead.
- `next()` advances an in-memory pointer through the flat byte array — **O(1), zero system calls** during iteration.

---

## Step 7 — K-mer Counting (`src/count.hpp` → `count_partition`)

For each superkmer delivered by the reader:

### 7a. Window initialisation — `init_packed_with_min(packed, min_pos)`

1. Decode the k packed bases into the `Kmer_Window` (one pass, O(k)).
2. Compute the ntHash of the minimizer l-mer by running `Roller<l>::init()` on the l bases starting at `min_pos` — **O(l) total**.
3. Store this hash as `precomp_nt_h_`; mark `use_precomp_ = true` so all subsequent `advance()` calls skip minimizer re-computation.

Net cost per superkmer: **O(k + l)** instead of the old O(k + l·sk_len) from a per-advance minimizer rescan.

### 7b. Prefetch

```cpp
table.prefetch(win);
```

Issues a CPU hardware prefetch for the primary hash bucket corresponding to the first k-mer. This hides the **~40 ns LLC miss** behind the O(l) minimizer computation above — the data arrives before the first `upsert`.

### 7c. First k-mer insert

```cpp
table.upsert(win, count++, token);
```

### 7d. Hot loop — remaining k-mers

```cpp
byte_ptr = packed + (k >> 2)
shift    = 6 − 2·(k & 3)           // initial bit offset, no division

for i = k … len−1:
    base  = (*byte_ptr >> shift) & 3  // 2-bit base, kache encoding
    shift -= 2
    if shift < 0: shift = 6; ++byte_ptr
    win.advance(base)
    table.upsert(win, count++, token)
```

- No ASCII conversion anywhere in the loop.
- Bit position tracked with `shift` counter — avoids division and modulo.
- `win.advance(base)` rolls the ntHash O(1) and updates the canonical k-mer state.

---

## Step 8 — Hash Table (`kache-hash`)

`Streaming_Kmer_Hash_Table<k, mt_=false, uint32_t, l>`:

- **Flat table** of fixed-size B-slot buckets plus an **overflow table** for high-load buckets.
- `mt_ = false` — no read-write locking; each partition's table is owned by exactly one thread.
- **`upsert(kmer, updater, default, token)`**: if the k-mer exists → apply `updater(old_count)` (i.e., `count + 1`); if new → insert `default` (1). Token is a per-thread reuse hint.
- **Resize** is triggered when: (a) main table load ≥ resize threshold, or (b) overflow insertion count ≥ overflow resize threshold. On resize: capacity doubles, flat + overflow tables are rehashed; the calling thread blocks; other partitions are unaffected.
- **`for_each(F)`**: iterates all (kmer, count) pairs in flat + overflow tables. Used by `write_counts`.

---

## Step 9 — Output (`src/count.hpp` → `write_counts`)

```
for each (kmer, count) in table.for_each():
    if count < ci or count > cx:  skip          // min/max abundance filters
    kmer.get_label() → ASCII string
    append "ACGTACGT…\tN\n" to thread-local chunk
    if chunk ≥ 1 MB:
        lock out_mutex → write to file → unlock
final flush of remaining chunk
```

- All Phase 2 threads write to the same output file behind a single mutex, batching in 1 MB chunks to amortise lock contention.
- Output order is not deterministic (hash table iteration order).

