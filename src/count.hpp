#pragma once

// Phase 2+3 — partition files → counted k-mers → output.
//
// Three independently modifiable bricks:
//
//   count_partition<k,l>   — read one partition file, fill a hash table.
//                            Swap to change the counting structure.
//
//   write_counts<k,l>      — drain one table to the output stream.
//                            Swap to change the output format.
//
//   count_and_write<k,l>   — parallel harness that calls both per partition.
//                            Swap to change the threading / scheduling strategy.

#include "Config.hpp"
#include "superkmer_io.hpp"

#include <fstream>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <utility>
#include <algorithm>

#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"


// ─── Counting brick ───────────────────────────────────────────────────────────
//
// Drain an already-open SuperkmerReader into the caller-provided hash table.
// The table must be private to the calling thread (mt_=false, no locking).
// Returns the total number of k-mer insertions (with multiplicity).

template <uint16_t k, uint16_t l>
uint64_t count_partition(
    SuperkmerReader&                                                      reader,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>::Token& token)
{
    uint64_t inserted = 0;
    while (reader.next()) {
        if (reader.size() < k) continue;

        kache_hash::Kmer_Window<k, l> win;
        win.init(reader.data());
        table.upsert(win, [](uint32_t v) { return v + 1; }, uint32_t(1), token);
        ++inserted;

        for (size_t i = k; i < reader.size(); ++i) {
            win.advance(reader[i]);
            table.upsert(win, [](uint32_t v) { return v + 1; }, uint32_t(1), token);
            ++inserted;
        }
    }
    return inserted;
}


// ─── Output brick ─────────────────────────────────────────────────────────────
//
// Iterate the table, apply ci/cx filters, and write "<kmer>\t<count>\n" lines
// in 1 MB batches under out_mutex.  chunk is a reusable per-thread buffer.
// Returns the number of k-mers written.
//
// Swap this function to change the output format (binary, CSV, different
// columns, post-processing, etc.) without touching counting logic.

template <uint16_t k, uint16_t l>
uint64_t write_counts(
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>& table,
    const Config&   cfg,
    std::string&    chunk,
    std::ofstream&  out,
    std::mutex&     out_mutex)
{
    constexpr size_t WRITE_BATCH = 1u << 20; // 1 MB

    const auto flush_chunk = [&]() {
        std::lock_guard<std::mutex> g(out_mutex);
        out.write(chunk.data(), chunk.size());
        chunk.clear();
    };

    std::string label;
    uint64_t written = 0;

    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        if (cnt < cfg.ci || cnt > cfg.cx) return;
        entry.first.get_label(label);
        chunk += label;
        chunk += '\t';
        chunk += std::to_string(cnt);
        chunk += '\n';
        ++written;
        if (chunk.size() >= WRITE_BATCH) flush_chunk();
    });
    if (!chunk.empty()) flush_chunk();
    return written;
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// Spawns min(num_threads, num_partitions) threads.  Each thread processes its
// assigned partitions round-robin: count_partition fills a private table, then
// write_counts drains it.  Returns {total_inserted, total_written}.

template <uint16_t k, uint16_t l>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    std::ofstream& out)
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};

    auto worker = [&](size_t tid) {
        typename table_t::Token token;
        std::string chunk;

        for (size_t p = tid; p < n_parts; p += n_threads) {
            // Open the partition file via mmap (SuperkmerReader maps on construction).
            // Use reader.file_size() to pre-size the hash table — avoids a separate
            // stat() call and reuses the already-mapped metadata.
            // Each byte encodes roughly one nucleotide; file_bytes / 4 is a
            // conservative upper bound on unique k-mers for this partition.
            const std::string part_path = cfg.work_dir + cfg.partition_prefix()
                                          + "_" + std::to_string(p) + ".superkmers";
            SuperkmerReader reader(part_path);
            const size_t init_sz = std::max(size_t(1) << 20, reader.file_size() / 4);
            table_t table(init_sz, 1); // private table, no locking

            const uint64_t ins = count_partition<k, l>(reader, table, token);
            total_inserted.fetch_add(ins, std::memory_order_relaxed);

            const uint64_t wrt = write_counts<k, l>(table, cfg, chunk, out, out_mutex);
            total_written.fetch_add(wrt, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    return { total_inserted.load(), total_written.load() };
}
