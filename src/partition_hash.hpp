#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// Partitioning: minimizer_hash % num_partitions.
//
// Parsing is delegated to SeqSource / GzInput + helicase SIMD parsers, which deliver
// ACTG-only chunks regardless of file format or compression.
//
// Internals:
//   extract_superkmers_from_actg<k, m, PartitionFn>  — pure ACTG sequence logic,
//                                                       no I/O, no threads.
//   partition_kmers_impl<k, m, PartitionFn>           — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "minimizer_window.hpp"
#include "seq_source.hpp"

#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <exception>
#include <memory>
#include <string_view>
#include <utility>
#include <array>
#include <cstring>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>


// Per-writer flush threshold: max(4 KB, budget_per_thread / n_parts).
// Disk mode: budget is RAM-proportional — larger buffers yield fewer, bigger writes.
// In-memory mode: fixed 64 MB budget (writes go to RAM, no I/O benefit).
inline size_t writer_flush_threshold(size_t n_parts, size_t budget_per_thread)
{
    constexpr size_t MIN_FLUSH_BYTES = 4u << 10;  // 4 KB minimum
    return std::max(MIN_FLUSH_BYTES, budget_per_thread / n_parts);
}

inline bool phase1_queue_stats_enabled()
{
    static const bool enabled = [] {
        const char* v = std::getenv("TUNA_PHASE1_QUEUE_STATS");
        return v && v[0] != '\0' && v[0] != '0';
    }();
    return enabled;
}

inline uint64_t phase1_now_ns()
{
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
}

inline double phase1_ns_to_s(uint64_t ns)
{
    return static_cast<double>(ns) / 1.0e9;
}

inline size_t phase1_env_size(const char* name, size_t fallback)
{
    const char* v = std::getenv(name);
    if (!v || v[0] == '\0') return fallback;
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(v, &end, 10);
    if (end == v || *end != '\0' || parsed == 0) return fallback;
    return static_cast<size_t>(parsed);
}

inline size_t phase1_gz_producer_threads(size_t n_files, size_t n_threads)
{
    if (n_threads <= 1) return 1;
    const size_t default_producers = std::max<size_t>(1, n_threads / 2);
    const size_t requested = phase1_env_size("TUNA_GZ_PRODUCERS", default_producers);
    return std::min(n_files, std::min(requested, n_threads - 1));
}

struct Phase1QueueThreadStats {
    uint64_t batches = 0;
    uint64_t records = 0;
    uint64_t bases = 0;
    uint64_t kmers = 0;
    uint64_t superkmers = 0;
    uint64_t wait_ns = 0;
};

inline void phase1_print_thread_stats(
    const char* label,
    const std::vector<Phase1QueueThreadStats>& stats)
{
    uint64_t batches = 0, records = 0, bases = 0, kmers = 0, superkmers = 0, wait_ns = 0;
    for (const auto& s : stats) {
        batches += s.batches;
        records += s.records;
        bases += s.bases;
        kmers += s.kmers;
        superkmers += s.superkmers;
        wait_ns += s.wait_ns;
    }
    std::cerr << "[phase1-queue] " << label
              << " total: batches=" << batches
              << " records=" << records
              << " bases=" << bases
              << " kmers=" << kmers
              << " superkmers=" << superkmers
              << " wait_s=" << phase1_ns_to_s(wait_ns) << "\n";
    for (size_t i = 0; i < stats.size(); ++i) {
        const auto& s = stats[i];
        std::cerr << "[phase1-queue] " << label << "[" << i << "]"
                  << " batches=" << s.batches
                  << " records=" << s.records
                  << " bases=" << s.bases
                  << " kmers=" << s.kmers
                  << " superkmers=" << s.superkmers
                  << " wait_s=" << phase1_ns_to_s(s.wait_ns) << "\n";
    }
}

struct Phase1AdaptiveWorkerStats {
    uint64_t decompress_tasks = 0;
    uint64_t parse_tasks = 0;
    uint64_t decompressed_chunks = 0;
    uint64_t decompressed_bytes = 0;
    uint64_t parsed_chunks = 0;
    uint64_t records = 0;
    uint64_t raw_bytes = 0;
    uint64_t bases = 0;
    uint64_t kmers = 0;
    uint64_t superkmers = 0;
    uint64_t wait_ns = 0;
};

inline void phase1_print_adaptive_stats(
    const char* label,
    const std::vector<Phase1AdaptiveWorkerStats>& stats,
    size_t max_queue_depth,
    size_t max_pending_raw_bytes,
    size_t low_raw_bytes,
    size_t high_raw_bytes,
    size_t opened_streams)
{
    uint64_t decompress_tasks = 0, parse_tasks = 0, decompressed_chunks = 0;
    uint64_t decompressed_bytes = 0, parsed_chunks = 0, records = 0, raw_bytes = 0;
    uint64_t bases = 0, kmers = 0, superkmers = 0, wait_ns = 0;
    for (const auto& s : stats) {
        decompress_tasks += s.decompress_tasks;
        parse_tasks += s.parse_tasks;
        decompressed_chunks += s.decompressed_chunks;
        decompressed_bytes += s.decompressed_bytes;
        parsed_chunks += s.parsed_chunks;
        records += s.records;
        raw_bytes += s.raw_bytes;
        bases += s.bases;
        kmers += s.kmers;
        superkmers += s.superkmers;
        wait_ns += s.wait_ns;
    }
    std::cerr << "[phase1-queue] " << label
              << " total: decompress_tasks=" << decompress_tasks
              << " parse_tasks=" << parse_tasks
              << " decompressed_chunks=" << decompressed_chunks
              << " decompressed_bytes=" << decompressed_bytes
              << " parsed_chunks=" << parsed_chunks
              << " records=" << records
              << " raw_bytes=" << raw_bytes
              << " bases=" << bases
              << " kmers=" << kmers
              << " superkmers=" << superkmers
              << " wait_s=" << phase1_ns_to_s(wait_ns)
              << " max_queue_depth=" << max_queue_depth
              << " max_pending_raw_bytes=" << max_pending_raw_bytes
              << " low_raw_bytes=" << low_raw_bytes
              << " high_raw_bytes=" << high_raw_bytes
              << " opened_streams=" << opened_streams << "\n";
    for (size_t i = 0; i < stats.size(); ++i) {
        const auto& s = stats[i];
        std::cerr << "[phase1-queue] " << label << "_worker[" << i << "]"
                  << " decompress_tasks=" << s.decompress_tasks
                  << " parse_tasks=" << s.parse_tasks
                  << " decompressed_chunks=" << s.decompressed_chunks
                  << " decompressed_bytes=" << s.decompressed_bytes
                  << " parsed_chunks=" << s.parsed_chunks
                  << " records=" << s.records
                  << " raw_bytes=" << s.raw_bytes
                  << " bases=" << s.bases
                  << " kmers=" << s.kmers
                  << " superkmers=" << s.superkmers
                  << " wait_s=" << phase1_ns_to_s(s.wait_ns) << "\n";
    }
}

// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only sequence, detect superkmer boundaries on minimizer hash
// changes, and append each superkmer to the corresponding writer.
// Splits on hash change (not partition change) so every superkmer has one minimizer.
// min_pos comes from MinimizerWindow::min_lmer_pos() — no extra scan needed.
//
// FlushFn: void(std::vector<SuperkmerWriter<k,m>>&, size_t partition_id)
// Called after each append; O(1) per superkmer.
template <uint16_t k, uint16_t m>
static constexpr uint16_t phase1_signature_len_v =
#ifdef TUNA_PHASE1_SIG_LEN
    (TUNA_PHASE1_SIG_LEN < k ? TUNA_PHASE1_SIG_LEN : m)
#else
    m
#endif
;

template <uint16_t k, uint16_t m, uint16_t sig_m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg_sig(
    const char* const                    seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, sig_m, (sig_m == m)>& min_it, // phase-1 signature iterator
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn,
    std::vector<uint8_t>&                kache_buf)   // kache-encoded sequence buffer (A=0,C=1,G=2,T=3)
{
    using hdr_t = sk_hdr_t<k, m>;  // superkmer header type (local alias)
    static_assert(sig_m < k, "phase-1 signature length must be strictly less than k");
    static_assert(sig_m <= 32, "phase-1 signature length must be <= 32");
    // Max value of hdr_t: flush guard prevents sk_len from ever overflowing hdr_t.
    // Rare long same-minimizer runs can exceed 2k-m, so use the serialized
    // header limit rather than the natural two-window span.
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    if (seq_len < k) return;

    // Pre-encode ASCII→kache once; bases are read multiple times across superkmers.
    kache_buf.resize(seq_len);
    for (size_t i = 0; i < seq_len; ++i)
        kache_buf[i] = ((uint8_t(seq[i]) >> 2) ^ (uint8_t(seq[i]) >> 1)) & 3u;

    min_it.reset(seq);   // initialises from ASCII (called once per sequence — cheap)
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = 0;
    if constexpr (sig_m == m)
        prev_min_pos = min_it.min_lmer_pos(); // absolute pos within seq
    size_t   pid          = partition_fn(prev_hash);  // partition id
    size_t   sk_start     = 0;                        // superkmer start position in current sequence

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance_kache(kache_buf[pos]);
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            hdr_t min_pos = NO_MIN;
            if constexpr (sig_m == m)
                min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
            writers[pid].append_kache(kache_buf.data() + sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            if constexpr (sig_m == m)
                prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    hdr_t min_pos = NO_MIN;
    if constexpr (sig_m == m)
        min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
    writers[pid].append_kache(kache_buf.data() + sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}

template <uint16_t k, uint16_t m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg(
    const char* const                    seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, m>&               min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn,
    std::vector<uint8_t>&                kache_buf)
{
    static constexpr uint16_t sig_m = phase1_signature_len_v<k, m>;
    if constexpr (sig_m == m) {
        extract_superkmers_from_actg_sig<k, m, m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), min_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn), kache_buf);
    } else {
        MinimizerWindow<k, sig_m, false> sig_it;
        extract_superkmers_from_actg_sig<k, m, sig_m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), sig_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn), kache_buf);
    }
}

template <uint16_t k, uint16_t m, uint16_t sig_m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_kache_sig(
    const uint8_t* const                 seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, sig_m, (sig_m == m)>& min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    using hdr_t = sk_hdr_t<k, m>;
    static_assert(sig_m < k, "phase-1 signature length must be strictly less than k");
    static_assert(sig_m <= 32, "phase-1 signature length must be <= 32");
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    if (seq_len < k) return;

    auto get_kache = [&](uint16_t i) noexcept -> uint8_t { return seq[i]; };

    min_it.reset_kache(get_kache);
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = 0;
    if constexpr (sig_m == m)
        prev_min_pos = min_it.min_lmer_pos();
    size_t   pid          = partition_fn(prev_hash);
    size_t   sk_start     = 0;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance_kache(seq[pos]);
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            hdr_t min_pos = NO_MIN;
            if constexpr (sig_m == m)
                min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
            writers[pid].append_kache(seq + sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            if constexpr (sig_m == m)
                prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    hdr_t min_pos = NO_MIN;
    if constexpr (sig_m == m)
        min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
    writers[pid].append_kache(seq + sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}

template <uint16_t k, uint16_t m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_kache(
    const uint8_t* const                 seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, m>&               min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    static constexpr uint16_t sig_m = phase1_signature_len_v<k, m>;
    if constexpr (sig_m == m) {
        extract_superkmers_from_kache_sig<k, m, m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), min_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    } else {
        MinimizerWindow<k, sig_m, false> sig_it;
        extract_superkmers_from_kache_sig<k, m, sig_m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), sig_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    }
}

template <uint16_t k, uint16_t m, uint16_t sig_m, typename PackedSeq, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_packed_nt_sig(
    const PackedSeq&                     seq,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, sig_m, (sig_m == m)>& min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    using hdr_t = sk_hdr_t<k, m>;
    static_assert(sig_m < k, "phase-1 signature length must be strictly less than k");
    static_assert(sig_m <= 32, "phase-1 signature length must be <= 32");
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    const size_t seq_len = seq.len();
    if (seq_len < k) return;

    auto get_nt = [&](size_t i) noexcept -> uint8_t { return seq.get(i); };

    min_it.reset_nt(get_nt);
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = 0;
    if constexpr (sig_m == m)
        prev_min_pos = min_it.min_lmer_pos();
    size_t   pid          = partition_fn(prev_hash);
    size_t   sk_start     = 0;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance_nt(get_nt(pos));
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            hdr_t min_pos = NO_MIN;
            if constexpr (sig_m == m)
                min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
            writers[pid].append_nt(get_nt, sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            if constexpr (sig_m == m)
                prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    hdr_t min_pos = NO_MIN;
    if constexpr (sig_m == m)
        min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
    writers[pid].append_nt(get_nt, sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}

template <uint16_t k, uint16_t m, typename PackedSeq, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_packed_nt(
    const PackedSeq&                     seq,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, m>&               min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    static constexpr uint16_t sig_m = phase1_signature_len_v<k, m>;
    if constexpr (sig_m == m) {
        extract_superkmers_from_packed_nt_sig<k, m, m>(
            seq, std::forward<PartitionFn>(partition_fn), min_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    } else {
        MinimizerWindow<k, sig_m, false> sig_it;
        extract_superkmers_from_packed_nt_sig<k, m, sig_m>(
            seq, std::forward<PartitionFn>(partition_fn), sig_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    }
}


using PackedDNAWord = __uint128_t;

struct PackedReadRecord {
    size_t        word_offset = 0;
    size_t        word_count  = 0;
    PackedDNAWord tail        = 0;
    size_t        len         = 0;
};

struct PackedReadBatch {
    std::vector<PackedDNAWord> words;
    std::vector<PackedReadRecord> records;
    size_t bases = 0;

    bool empty() const noexcept { return records.empty(); }

    void reserve(size_t max_records, size_t target_bases)
    {
        records.reserve(max_records);
        words.reserve((target_bases + 63u) / 64u + max_records);
    }

    void append(const helicase::PackedDNA& dna)
    {
        auto [dna_words, tail] = dna.bits();
        PackedReadRecord rec;
        rec.word_offset = words.size();
        rec.word_count  = dna_words.size();
        rec.tail        = tail;
        rec.len         = dna.len();
        words.insert(words.end(), dna_words.begin(), dna_words.end());
        records.push_back(rec);
        bases += rec.len;
    }
};

struct PackedReadView {
    const PackedDNAWord* words = nullptr;
    size_t               word_count = 0;
    PackedDNAWord        tail = 0;
    size_t               len_bases = 0;

    size_t len() const noexcept { return len_bases; }

    uint8_t get(size_t i) const noexcept
    {
        const size_t word = i / 64u;
        const size_t bit  = i % 64u;
        const PackedDNAWord value = word < word_count ? words[word] : tail;
        return static_cast<uint8_t>((value >> (2u * bit)) & 0b11u);
    }
};

inline constexpr auto make_nt_byte_to_kache_lut()
{
    std::array<std::array<uint8_t, 4>, 256> lut{};
    for (size_t x = 0; x < lut.size(); ++x) {
        for (size_t i = 0; i < 4; ++i) {
            const uint8_t nt = static_cast<uint8_t>((x >> (2u * i)) & 0b11u);
            lut[x][i] = static_cast<uint8_t>(nt ^ (nt >> 1));
        }
    }
    return lut;
}

inline void decode_packed_nt_to_kache(const PackedReadView& rec, std::vector<uint8_t>& out)
{
    static constexpr auto LUT = make_nt_byte_to_kache_lut();

    out.resize(rec.len());
    uint8_t* dst = out.data();
    size_t remaining = rec.len();

    auto decode_word = [&](PackedDNAWord word) {
        constexpr size_t BASES_PER_WORD = 64;
        constexpr size_t BASES_PER_BYTE = 4;
        const size_t n = std::min(remaining, BASES_PER_WORD);
        const size_t full_groups = n / BASES_PER_BYTE;
        for (size_t g = 0; g < full_groups; ++g) {
            const auto& bases = LUT[static_cast<uint8_t>(word >> (8u * g))];
            std::memcpy(dst, bases.data(), BASES_PER_BYTE);
            dst += BASES_PER_BYTE;
        }
        for (size_t i = full_groups * BASES_PER_BYTE; i < n; ++i) {
            const uint8_t nt = static_cast<uint8_t>((word >> (2u * i)) & 0b11u);
            *dst++ = static_cast<uint8_t>(nt ^ (nt >> 1));
        }
        remaining -= n;
    };

    for (size_t w = 0; w < rec.word_count && remaining > 0; ++w)
        decode_word(rec.words[w]);
    if (remaining > 0)
        decode_word(rec.tail);
}

inline void decode_packed_nt_to_kache(const helicase::PackedDNA& dna, std::vector<uint8_t>& out)
{
    auto [words, tail] = dna.bits();
    PackedReadView rec{words.data(), words.size(), tail, dna.len()};
    decode_packed_nt_to_kache(rec, out);
}

struct RawFastqChunk {
    std::vector<uint8_t> data;

    bool empty() const noexcept { return data.empty(); }
};

inline uint8_t gz_first_byte(const std::string& path)
{
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) throw std::runtime_error("Cannot open: " + path);
    gzbuffer(gz, 1u << 20);
    uint8_t c = 0;
    const int n = gzread(gz, &c, 1);
    gzclose(gz);
    if (n != 1) throw std::runtime_error("Cannot read: " + path);
    return c;
}

inline bool phase1_adaptive_all_gz_fastq(const Config& cfg, bool all_gz)
{
    if (!cfg.phase1_adaptive || !all_gz) return false;
    for (const auto& f : cfg.input_files) {
        if (gz_first_byte(f) != '@')
            return false;
    }
    return true;
}

class RawGzFastqChunker {
    static constexpr size_t READ_BYTES = 4u << 20;

    gzFile gz_ = nullptr;
    std::vector<uint8_t> buf_;
    size_t target_bytes_ = 0;
    bool eof_ = false;

    static size_t last_record_boundary(const std::vector<uint8_t>& data) noexcept
    {
        size_t line_count = 0;
        size_t boundary = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                ++line_count;
                if ((line_count & 3u) == 0)
                    boundary = i + 1;
            }
        }
        return boundary;
    }

    void read_more()
    {
        const size_t old_size = buf_.size();
        buf_.resize(old_size + READ_BYTES);
        const int n = gzread(gz_, buf_.data() + old_size, READ_BYTES);
        if (n < 0) {
            int errnum = 0;
            const char* msg = gzerror(gz_, &errnum);
            throw std::runtime_error(std::string("gzip read failed: ") + (msg ? msg : "unknown error"));
        }
        if (n == 0) {
            eof_ = true;
            buf_.resize(old_size);
        } else {
            buf_.resize(old_size + static_cast<size_t>(n));
        }
    }

public:
    RawGzFastqChunker(const std::string& path, size_t target_bytes)
        : target_bytes_(target_bytes)
    {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) throw std::runtime_error("Cannot open: " + path);
        gzbuffer(gz_, 1u << 20);
        buf_.reserve(target_bytes_ + READ_BYTES);
    }

    ~RawGzFastqChunker()
    {
        if (gz_) gzclose(gz_);
    }

    RawGzFastqChunker(const RawGzFastqChunker&) = delete;
    RawGzFastqChunker& operator=(const RawGzFastqChunker&) = delete;

    bool next(RawFastqChunk& out)
    {
        out.data.clear();

        while (!eof_ && buf_.size() < target_bytes_)
            read_more();

        if (buf_.empty())
            return false;

        size_t boundary = last_record_boundary(buf_);
        while (!eof_ && boundary == 0) {
            read_more();
            boundary = last_record_boundary(buf_);
        }

        if (eof_) {
            out.data = std::move(buf_);
            buf_.clear();
            return !out.data.empty();
        }

        RawFastqChunk chunk;
        chunk.data = std::move(buf_);
        std::vector<uint8_t> tail;
        if (boundary < chunk.data.size())
            tail.assign(chunk.data.begin() + static_cast<std::ptrdiff_t>(boundary), chunk.data.end());
        chunk.data.resize(boundary);
        buf_ = std::move(tail);
        out = std::move(chunk);
        return !out.data.empty();
    }
};

class AsyncPartitionWriters {
    struct Shard {
        std::mutex mtx;
        std::condition_variable cv;
        std::deque<SuperkmerWriteBlock> queue;
        size_t queued_bytes = 0;
        bool done = false;
        uint64_t enqueued_blocks = 0;
        uint64_t enqueued_bytes = 0;
        uint64_t written_blocks = 0;
        uint64_t written_bytes = 0;
        uint64_t enqueue_wait_ns = 0;
        uint64_t worker_wait_ns = 0;
        size_t max_queued_bytes = 0;
        size_t max_queue_depth = 0;
    };

    std::vector<SuperkmerBucketFile>& buckets_;
    std::vector<std::unique_ptr<Shard>> shards_;
    std::vector<std::thread> threads_;
#ifdef TUNA_LZ4_BUCKETS
    std::vector<uint8_t> initialized_;
    bool lz4_buckets_ = false;
#endif
    size_t max_queue_bytes_;
    bool stats_enabled_ = false;
    std::exception_ptr error_ = nullptr;
    std::mutex error_mutex_;

    size_t shard_for(size_t partition) const noexcept
    {
        return partition % shards_.size();
    }

    void record_error()
    {
        std::lock_guard<std::mutex> lk(error_mutex_);
        if (!error_) error_ = std::current_exception();
    }

    void worker(size_t shard_id)
    {
        auto& shard = *shards_[shard_id];
        bool failed = false;
        while (true) {
            SuperkmerWriteBlock block;
            {
                std::unique_lock<std::mutex> lk(shard.mtx);
                const uint64_t wait_start = stats_enabled_ ? phase1_now_ns() : 0;
                shard.cv.wait(lk, [&]{ return shard.done || !shard.queue.empty(); });
                if (stats_enabled_) shard.worker_wait_ns += phase1_now_ns() - wait_start;
                if (shard.queue.empty()) {
                    if (shard.done) break;
                    continue;
                }
                block = std::move(shard.queue.front());
                shard.queue.pop_front();
                shard.queued_bytes -= block.size;
                if (stats_enabled_) {
                    ++shard.written_blocks;
                    shard.written_bytes += block.size;
                }
            }
            shard.cv.notify_all();

            if (failed) continue;
            try {
                auto& out = buckets_[block.partition];
#ifdef TUNA_LZ4_BUCKETS
                if (lz4_buckets_ && !initialized_[block.partition]) {
                    tuna_superkmer_detail::write_lz4_bucket_magic(out);
                    initialized_[block.partition] = 1;
                }
#endif
                tuna_superkmer_detail::write_superkmer_bucket_block(
                    out, block.data, block.size,
#ifdef TUNA_LZ4_BUCKETS
                    lz4_buckets_
#else
                    false
#endif
                );
                if (!out) throw std::runtime_error("tuna: failed while writing partition file");
            } catch (...) {
                failed = true;
                record_error();
            }
        }
    }

public:
    AsyncPartitionWriters(std::vector<SuperkmerBucketFile>& buckets,
                          size_t n_shards,
                          size_t max_queue_bytes,
                          bool lz4_buckets = false)
        : buckets_(buckets),
#ifdef TUNA_LZ4_BUCKETS
          initialized_(buckets_.size(), 0),
          lz4_buckets_(lz4_buckets),
#endif
          max_queue_bytes_(max_queue_bytes)
    {
        stats_enabled_ = phase1_queue_stats_enabled();
        n_shards = std::max<size_t>(1, std::min(n_shards, buckets_.size()));
        shards_.reserve(n_shards);
        for (size_t i = 0; i < n_shards; ++i)
            shards_.push_back(std::make_unique<Shard>());
        threads_.reserve(n_shards);
        for (size_t i = 0; i < n_shards; ++i)
            threads_.emplace_back([this, i]{ worker(i); });
    }

    AsyncPartitionWriters(const AsyncPartitionWriters&) = delete;
    AsyncPartitionWriters& operator=(const AsyncPartitionWriters&) = delete;

    void enqueue(SuperkmerWriteBlock&& block)
    {
        if (block.size == 0) return;
        auto& shard = *shards_[shard_for(block.partition)];
        {
            std::unique_lock<std::mutex> lk(shard.mtx);
            const uint64_t wait_start = stats_enabled_ ? phase1_now_ns() : 0;
            shard.cv.wait(lk, [&]{
                return shard.done ||
                       shard.queued_bytes + block.size <= max_queue_bytes_;
            });
            if (stats_enabled_) shard.enqueue_wait_ns += phase1_now_ns() - wait_start;
            if (shard.done) return;
            shard.queued_bytes += block.size;
            if (stats_enabled_) {
                ++shard.enqueued_blocks;
                shard.enqueued_bytes += block.size;
                shard.max_queued_bytes = std::max(shard.max_queued_bytes, shard.queued_bytes);
                shard.max_queue_depth = std::max(shard.max_queue_depth, shard.queue.size() + 1);
            }
            shard.queue.emplace_back(std::move(block));
        }
        shard.cv.notify_one();
    }

    void finish()
    {
        for (auto& shard_ptr : shards_) {
            auto& shard = *shard_ptr;
            {
                std::lock_guard<std::mutex> lk(shard.mtx);
                shard.done = true;
            }
            shard.cv.notify_all();
        }
        for (auto& th : threads_)
            if (th.joinable()) th.join();
        if (stats_enabled_) {
            uint64_t enq_blocks = 0, enq_bytes = 0, wr_blocks = 0, wr_bytes = 0;
            uint64_t enq_wait_ns = 0, worker_wait_ns = 0;
            size_t max_bytes = 0, max_depth = 0;
            for (size_t i = 0; i < shards_.size(); ++i) {
                auto& s = *shards_[i];
                enq_blocks += s.enqueued_blocks;
                enq_bytes += s.enqueued_bytes;
                wr_blocks += s.written_blocks;
                wr_bytes += s.written_bytes;
                enq_wait_ns += s.enqueue_wait_ns;
                worker_wait_ns += s.worker_wait_ns;
                max_bytes = std::max(max_bytes, s.max_queued_bytes);
                max_depth = std::max(max_depth, s.max_queue_depth);
                std::cerr << "[phase1-queue] writer_shard[" << i << "]"
                          << " enqueued_blocks=" << s.enqueued_blocks
                          << " enqueued_bytes=" << s.enqueued_bytes
                          << " written_blocks=" << s.written_blocks
                          << " written_bytes=" << s.written_bytes
                          << " enqueue_wait_s=" << phase1_ns_to_s(s.enqueue_wait_ns)
                          << " worker_wait_s=" << phase1_ns_to_s(s.worker_wait_ns)
                          << " max_queued_bytes=" << s.max_queued_bytes
                          << " max_queue_depth=" << s.max_queue_depth << "\n";
            }
            std::cerr << "[phase1-queue] writer total:"
                      << " shards=" << shards_.size()
                      << " enqueued_blocks=" << enq_blocks
                      << " enqueued_bytes=" << enq_bytes
                      << " written_blocks=" << wr_blocks
                      << " written_bytes=" << wr_bytes
                      << " enqueue_wait_s=" << phase1_ns_to_s(enq_wait_ns)
                      << " worker_wait_s=" << phase1_ns_to_s(worker_wait_ns)
                      << " max_queued_bytes=" << max_bytes
                      << " max_queue_depth=" << max_depth << "\n";
        }
        if (error_) std::rethrow_exception(error_);
    }
};


struct AdaptiveRawChunk {
    RawFastqChunk chunk;
};

struct AdaptiveGzStream {
    std::unique_ptr<RawGzFastqChunker> chunker;
    size_t file_idx = 0;
};

using AdaptivePackedFastqParser = helicase::FastqParser<HELICASE_ACTG_PACKED, GzInput>;

struct AdaptivePackedGzStream {
    std::unique_ptr<AdaptivePackedFastqParser> parser;
    size_t file_idx = 0;
};

// Adaptive gz FASTQ path: all workers share a pool of gzip streams and a raw
// FASTQ chunk queue.  Workers choose decompression or parse/partition work from
// queue pressure instead of a fixed producer/consumer split.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_adaptive_gz_fastq_pc(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
#ifdef TUNA_FASTQ_CHUNK_MB
    constexpr size_t TARGET_CHUNK_BYTES = static_cast<size_t>(TUNA_FASTQ_CHUNK_MB) << 20;
#else
    constexpr size_t TARGET_CHUNK_BYTES = 4u << 20;
#endif

    const size_t n_threads = std::max<size_t>(1, static_cast<size_t>(cfg.num_threads));
    const size_t n_files = cfg.input_files.size();
    const size_t n_parts = cfg.num_partitions;
    const size_t initial_streams = std::min(n_files, n_threads);
    const size_t low_raw_bytes = TARGET_CHUNK_BYTES * std::max<size_t>(2, n_threads / 2);
    const size_t high_raw_bytes = TARGET_CHUNK_BYTES * std::max<size_t>(4, n_threads * 2);

    std::vector<AdaptiveGzStream> streams(initial_streams);
    std::deque<size_t> ready_streams;
    for (size_t i = 0; i < initial_streams; ++i) {
        streams[i].file_idx = i;
        streams[i].chunker = std::make_unique<RawGzFastqChunker>(cfg.input_files[i], TARGET_CHUNK_BYTES);
        ready_streams.push_back(i);
    }
    size_t next_file = initial_streams;
    size_t active_streams = initial_streams;
    size_t opened_streams = initial_streams;

    std::deque<AdaptiveRawChunk> raw_chunks;
    size_t pending_raw_bytes = 0;
    size_t max_queue_depth = 0;
    size_t max_pending_raw_bytes = 0;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    std::vector<Phase1AdaptiveWorkerStats> worker_stats(n_threads);

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);

    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto record_error = [&]() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
    };

    auto worker = [&](size_t worker_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                RawFastqChunk batch;
                size_t stream_id = 0;
                bool do_parse = false;
                bool do_decompress = false;

                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return stop.load(std::memory_order_relaxed) ||
                               !raw_chunks.empty() ||
                               (!ready_streams.empty() && pending_raw_bytes < high_raw_bytes) ||
                               active_streams == 0;
                    });
                    if (queue_stats) worker_stats[worker_id].wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;

                    const bool can_parse = !raw_chunks.empty();
                    const bool can_decompress = !ready_streams.empty() && pending_raw_bytes < high_raw_bytes;
                    if (can_parse && (!can_decompress || pending_raw_bytes >= low_raw_bytes)) {
                        batch = std::move(raw_chunks.front().chunk);
                        raw_chunks.pop_front();
                        pending_raw_bytes -= batch.data.size();
                        do_parse = true;
                    } else if (can_decompress) {
                        stream_id = ready_streams.front();
                        ready_streams.pop_front();
                        do_decompress = true;
                    } else if (can_parse) {
                        batch = std::move(raw_chunks.front().chunk);
                        raw_chunks.pop_front();
                        pending_raw_bytes -= batch.data.size();
                        do_parse = true;
                    } else if (active_streams == 0) {
                        break;
                    } else {
                        continue;
                    }
                }
                q_cv.notify_all();

                if (do_decompress) {
                    ++worker_stats[worker_id].decompress_tasks;
                    RawFastqChunk produced;
                    const bool have_chunk = streams[stream_id].chunker->next(produced);
                    if (have_chunk) {
                        const size_t produced_bytes = produced.data.size();
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            raw_chunks.push_back(AdaptiveRawChunk{std::move(produced)});
                            pending_raw_bytes += produced_bytes;
                            ready_streams.push_back(stream_id);
                            max_queue_depth = std::max(max_queue_depth, raw_chunks.size());
                            max_pending_raw_bytes = std::max(max_pending_raw_bytes, pending_raw_bytes);
                        }
                        ++worker_stats[worker_id].decompressed_chunks;
                        worker_stats[worker_id].decompressed_bytes += produced_bytes;
                        q_cv.notify_all();
                    } else {
                        streams[stream_id].chunker.reset();
                        std::string next_path;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (next_file < n_files) {
                                next_path = cfg.input_files[next_file++];
                                ++opened_streams;
                            } else {
                                --active_streams;
                            }
                        }
                        if (!next_path.empty()) {
                            auto next_chunker = std::make_unique<RawGzFastqChunker>(next_path, TARGET_CHUNK_BYTES);
                            {
                                std::lock_guard<std::mutex> lk(q_mutex);
                                streams[stream_id].chunker = std::move(next_chunker);
                                ready_streams.push_back(stream_id);
                            }
                        }
                        q_cv.notify_all();
                    }
                    continue;
                }

                if (do_parse) {
                    ++worker_stats[worker_id].parse_tasks;
                    ++worker_stats[worker_id].parsed_chunks;
                    worker_stats[worker_id].raw_bytes += batch.data.size();
                    helicase::FastqParser<HELICASE_ACTG_PACKED, helicase::SliceInput> parser(
                        batch.data.data(), batch.data.size());
                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        const size_t len = parser.get_dna_len();
                        if (len < k) continue;
#ifdef TUNA_PACKED_NT_PHASE1
                        auto dna = parser.get_dna_packed();
                        auto [words, tail] = dna.bits();
                        PackedReadView rec{words.data(), words.size(), tail, dna.len()};
                        extract_superkmers_from_packed_nt<k, m>(
                            rec, partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#else
                        decode_packed_nt_to_kache(parser.get_dna_packed(), kache_buf);
                        extract_superkmers_from_kache<k, m>(
                            kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#endif
                        ++local_seqs;
                        ++worker_stats[worker_id].records;
                        worker_stats[worker_id].bases += len;
                    }
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            worker_stats[worker_id].kmers += local_kmers;
            worker_stats[worker_id].superkmers += local_superkmers;
        } catch (...) {
            record_error();
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (worker_error) std::rethrow_exception(worker_error);
    if (queue_stats) {
        phase1_print_adaptive_stats(
            "adaptive_gz_fastq", worker_stats, max_queue_depth,
            max_pending_raw_bytes, low_raw_bytes, high_raw_bytes, opened_streams);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_adaptive_gz_fastq_pc(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
#ifdef TUNA_FASTQ_CHUNK_MB
    constexpr size_t TARGET_CHUNK_BYTES = static_cast<size_t>(TUNA_FASTQ_CHUNK_MB) << 20;
#else
    constexpr size_t TARGET_CHUNK_BYTES = 4u << 20;
#endif

    const size_t n_threads = std::max<size_t>(1, static_cast<size_t>(cfg.num_threads));
    const size_t n_files = cfg.input_files.size();
    const size_t n_parts = cfg.num_partitions;
    const size_t initial_streams = std::min(n_files, n_threads);
    const size_t low_raw_bytes = TARGET_CHUNK_BYTES * std::max<size_t>(2, n_threads / 2);
    const size_t high_raw_bytes = TARGET_CHUNK_BYTES * std::max<size_t>(4, n_threads * 2);

    std::vector<AdaptiveGzStream> streams(initial_streams);
    std::deque<size_t> ready_streams;
    for (size_t i = 0; i < initial_streams; ++i) {
        streams[i].file_idx = i;
        streams[i].chunker = std::make_unique<RawGzFastqChunker>(cfg.input_files[i], TARGET_CHUNK_BYTES);
        ready_streams.push_back(i);
    }
    size_t next_file = initial_streams;
    size_t active_streams = initial_streams;
    size_t opened_streams = initial_streams;

    std::deque<AdaptiveRawChunk> raw_chunks;
    size_t pending_raw_bytes = 0;
    size_t max_queue_depth = 0;
    size_t max_pending_raw_bytes = 0;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    std::vector<Phase1AdaptiveWorkerStats> worker_stats(n_threads);
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto record_error = [&]() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
    };

    auto worker = [&](size_t worker_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                RawFastqChunk batch;
                size_t stream_id = 0;
                bool do_parse = false;
                bool do_decompress = false;

                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return stop.load(std::memory_order_relaxed) ||
                               !raw_chunks.empty() ||
                               (!ready_streams.empty() && pending_raw_bytes < high_raw_bytes) ||
                               active_streams == 0;
                    });
                    if (queue_stats) worker_stats[worker_id].wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;

                    const bool can_parse = !raw_chunks.empty();
                    const bool can_decompress = !ready_streams.empty() && pending_raw_bytes < high_raw_bytes;
                    if (can_parse && (!can_decompress || pending_raw_bytes >= low_raw_bytes)) {
                        batch = std::move(raw_chunks.front().chunk);
                        raw_chunks.pop_front();
                        pending_raw_bytes -= batch.data.size();
                        do_parse = true;
                    } else if (can_decompress) {
                        stream_id = ready_streams.front();
                        ready_streams.pop_front();
                        do_decompress = true;
                    } else if (can_parse) {
                        batch = std::move(raw_chunks.front().chunk);
                        raw_chunks.pop_front();
                        pending_raw_bytes -= batch.data.size();
                        do_parse = true;
                    } else if (active_streams == 0) {
                        break;
                    } else {
                        continue;
                    }
                }
                q_cv.notify_all();

                if (do_decompress) {
                    ++worker_stats[worker_id].decompress_tasks;
                    RawFastqChunk produced;
                    const bool have_chunk = streams[stream_id].chunker->next(produced);
                    if (have_chunk) {
                        const size_t produced_bytes = produced.data.size();
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            raw_chunks.push_back(AdaptiveRawChunk{std::move(produced)});
                            pending_raw_bytes += produced_bytes;
                            ready_streams.push_back(stream_id);
                            max_queue_depth = std::max(max_queue_depth, raw_chunks.size());
                            max_pending_raw_bytes = std::max(max_pending_raw_bytes, pending_raw_bytes);
                        }
                        ++worker_stats[worker_id].decompressed_chunks;
                        worker_stats[worker_id].decompressed_bytes += produced_bytes;
                        q_cv.notify_all();
                    } else {
                        streams[stream_id].chunker.reset();
                        std::string next_path;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (next_file < n_files) {
                                next_path = cfg.input_files[next_file++];
                                ++opened_streams;
                            } else {
                                --active_streams;
                            }
                        }
                        if (!next_path.empty()) {
                            auto next_chunker = std::make_unique<RawGzFastqChunker>(next_path, TARGET_CHUNK_BYTES);
                            {
                                std::lock_guard<std::mutex> lk(q_mutex);
                                streams[stream_id].chunker = std::move(next_chunker);
                                ready_streams.push_back(stream_id);
                            }
                        }
                        q_cv.notify_all();
                    }
                    continue;
                }

                if (do_parse) {
                    ++worker_stats[worker_id].parse_tasks;
                    ++worker_stats[worker_id].parsed_chunks;
                    worker_stats[worker_id].raw_bytes += batch.data.size();
                    helicase::FastqParser<HELICASE_ACTG_PACKED, helicase::SliceInput> parser(
                        batch.data.data(), batch.data.size());
                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        const size_t len = parser.get_dna_len();
                        if (len < k) continue;
#ifdef TUNA_PACKED_NT_PHASE1
                        auto dna = parser.get_dna_packed();
                        auto [words, tail] = dna.bits();
                        PackedReadView rec{words.data(), words.size(), tail, dna.len()};
                        extract_superkmers_from_packed_nt<k, m>(
                            rec, partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#else
                        decode_packed_nt_to_kache(parser.get_dna_packed(), kache_buf);
                        extract_superkmers_from_kache<k, m>(
                            kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#endif
                        ++local_seqs;
                        ++worker_stats[worker_id].records;
                        worker_stats[worker_id].bases += len;
                    }
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            worker_stats[worker_id].kmers += local_kmers;
            worker_stats[worker_id].superkmers += local_superkmers;
        } catch (...) {
            record_error();
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);
    if (queue_stats) {
        phase1_print_adaptive_stats(
            "mem_adaptive_gz_fastq", worker_stats, max_queue_depth,
            max_pending_raw_bytes, low_raw_bytes, high_raw_bytes, opened_streams);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_adaptive_packed_gz_fastq_pc(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_threads = std::max<size_t>(1, static_cast<size_t>(cfg.num_threads));
    const size_t n_files = cfg.input_files.size();
    const size_t n_parts = cfg.num_partitions;
    const size_t initial_streams = std::min(n_files, n_threads);
    const size_t low_batch_bases = TARGET_BATCH_BASES * std::max<size_t>(1, n_threads / 2);
    const size_t high_batch_bases = TARGET_BATCH_BASES * std::max<size_t>(2, n_threads);

    std::vector<AdaptivePackedGzStream> streams(initial_streams);
    std::deque<size_t> ready_streams;
    for (size_t i = 0; i < initial_streams; ++i) {
        streams[i].file_idx = i;
        GzInput inp(cfg.input_files[i]);
        streams[i].parser = std::make_unique<AdaptivePackedFastqParser>(std::move(inp));
        ready_streams.push_back(i);
    }
    size_t next_file = initial_streams;
    size_t active_streams = initial_streams;
    size_t opened_streams = initial_streams;

    std::deque<Batch> batches;
    size_t pending_bases = 0;
    size_t max_queue_depth = 0;
    size_t max_pending_bases = 0;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    std::vector<Phase1AdaptiveWorkerStats> worker_stats(n_threads);

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto record_error = [&]() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
    };

    auto worker = [&](size_t worker_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                Batch batch;
                size_t stream_id = 0;
                bool do_parse = false;
                bool do_produce = false;

                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return stop.load(std::memory_order_relaxed) ||
                               !batches.empty() ||
                               (!ready_streams.empty() && pending_bases < high_batch_bases) ||
                               active_streams == 0;
                    });
                    if (queue_stats) worker_stats[worker_id].wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;

                    const bool can_parse = !batches.empty();
                    const bool can_produce = !ready_streams.empty() && pending_bases < high_batch_bases;
                    if (can_parse && (!can_produce || pending_bases >= low_batch_bases)) {
                        batch = std::move(batches.front());
                        batches.pop_front();
                        pending_bases -= batch.bases;
                        do_parse = true;
                    } else if (can_produce) {
                        stream_id = ready_streams.front();
                        ready_streams.pop_front();
                        do_produce = true;
                    } else if (can_parse) {
                        batch = std::move(batches.front());
                        batches.pop_front();
                        pending_bases -= batch.bases;
                        do_parse = true;
                    } else if (active_streams == 0) {
                        break;
                    } else {
                        continue;
                    }
                }
                q_cv.notify_all();

                if (do_produce) {
                    ++worker_stats[worker_id].decompress_tasks;
                    Batch produced;
                    produced.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                    size_t produced_bases = 0;
                    bool eof = false;
                    auto& parser = *streams[stream_id].parser;
                    while (!stop.load(std::memory_order_relaxed) &&
                           produced_bases < TARGET_BATCH_BASES &&
                           produced.records.size() < MAX_BATCH_ITEMS) {
                        if (!parser.next()) {
                            eof = true;
                            break;
                        }
                        const size_t len = parser.get_dna_len();
                        if (len >= k) {
                            produced.append(parser.get_dna_packed());
                            produced_bases += len;
                        }
                    }

                    if (!produced.empty()) {
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            pending_bases += produced.bases;
                            batches.push_back(std::move(produced));
                            max_queue_depth = std::max(max_queue_depth, batches.size());
                            max_pending_bases = std::max(max_pending_bases, pending_bases);
                        }
                        ++worker_stats[worker_id].decompressed_chunks;
                        worker_stats[worker_id].decompressed_bytes += produced_bases;
                    }

                    if (!eof && !stop.load(std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        ready_streams.push_back(stream_id);
                    } else {
                        streams[stream_id].parser.reset();
                        std::string next_path;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (next_file < n_files) {
                                next_path = cfg.input_files[next_file++];
                                ++opened_streams;
                            } else {
                                --active_streams;
                            }
                        }
                        if (!next_path.empty()) {
                            GzInput inp(next_path);
                            auto next_parser = std::make_unique<AdaptivePackedFastqParser>(std::move(inp));
                            {
                                std::lock_guard<std::mutex> lk(q_mutex);
                                streams[stream_id].parser = std::move(next_parser);
                                ready_streams.push_back(stream_id);
                            }
                        }
                    }
                    q_cv.notify_all();
                    continue;
                }

                if (do_parse) {
                    ++worker_stats[worker_id].parse_tasks;
                    ++worker_stats[worker_id].parsed_chunks;
                    worker_stats[worker_id].records += batch.records.size();
                    worker_stats[worker_id].bases += batch.bases;
                    for (const auto& rec : batch.records) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        PackedReadView chunk{
                            batch.words.data() + rec.word_offset,
                            rec.word_count,
                            rec.tail,
                            rec.len
                        };
#ifdef TUNA_PACKED_NT_PHASE1
                        extract_superkmers_from_packed_nt<k, m>(
                            chunk, partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#else
                        decode_packed_nt_to_kache(chunk, kache_buf);
                        extract_superkmers_from_kache<k, m>(
                            kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#endif
                        ++local_seqs;
                    }
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            worker_stats[worker_id].kmers += local_kmers;
            worker_stats[worker_id].superkmers += local_superkmers;
        } catch (...) {
            record_error();
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (worker_error) std::rethrow_exception(worker_error);
    if (queue_stats) {
        phase1_print_adaptive_stats(
            "adaptive_packed_gz_fastq", worker_stats, max_queue_depth,
            max_pending_bases, low_batch_bases, high_batch_bases, opened_streams);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_adaptive_packed_gz_fastq_pc(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_threads = std::max<size_t>(1, static_cast<size_t>(cfg.num_threads));
    const size_t n_files = cfg.input_files.size();
    const size_t n_parts = cfg.num_partitions;
    const size_t initial_streams = std::min(n_files, n_threads);
    const size_t low_batch_bases = TARGET_BATCH_BASES * std::max<size_t>(1, n_threads / 2);
    const size_t high_batch_bases = TARGET_BATCH_BASES * std::max<size_t>(2, n_threads);

    std::vector<AdaptivePackedGzStream> streams(initial_streams);
    std::deque<size_t> ready_streams;
    for (size_t i = 0; i < initial_streams; ++i) {
        streams[i].file_idx = i;
        GzInput inp(cfg.input_files[i]);
        streams[i].parser = std::make_unique<AdaptivePackedFastqParser>(std::move(inp));
        ready_streams.push_back(i);
    }
    size_t next_file = initial_streams;
    size_t active_streams = initial_streams;
    size_t opened_streams = initial_streams;

    std::deque<Batch> batches;
    size_t pending_bases = 0;
    size_t max_queue_depth = 0;
    size_t max_pending_bases = 0;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    std::vector<Phase1AdaptiveWorkerStats> worker_stats(n_threads);
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto record_error = [&]() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
    };

    auto worker = [&](size_t worker_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                Batch batch;
                size_t stream_id = 0;
                bool do_parse = false;
                bool do_produce = false;

                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return stop.load(std::memory_order_relaxed) ||
                               !batches.empty() ||
                               (!ready_streams.empty() && pending_bases < high_batch_bases) ||
                               active_streams == 0;
                    });
                    if (queue_stats) worker_stats[worker_id].wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;

                    const bool can_parse = !batches.empty();
                    const bool can_produce = !ready_streams.empty() && pending_bases < high_batch_bases;
                    if (can_parse && (!can_produce || pending_bases >= low_batch_bases)) {
                        batch = std::move(batches.front());
                        batches.pop_front();
                        pending_bases -= batch.bases;
                        do_parse = true;
                    } else if (can_produce) {
                        stream_id = ready_streams.front();
                        ready_streams.pop_front();
                        do_produce = true;
                    } else if (can_parse) {
                        batch = std::move(batches.front());
                        batches.pop_front();
                        pending_bases -= batch.bases;
                        do_parse = true;
                    } else if (active_streams == 0) {
                        break;
                    } else {
                        continue;
                    }
                }
                q_cv.notify_all();

                if (do_produce) {
                    ++worker_stats[worker_id].decompress_tasks;
                    Batch produced;
                    produced.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                    size_t produced_bases = 0;
                    bool eof = false;
                    auto& parser = *streams[stream_id].parser;
                    while (!stop.load(std::memory_order_relaxed) &&
                           produced_bases < TARGET_BATCH_BASES &&
                           produced.records.size() < MAX_BATCH_ITEMS) {
                        if (!parser.next()) {
                            eof = true;
                            break;
                        }
                        const size_t len = parser.get_dna_len();
                        if (len >= k) {
                            produced.append(parser.get_dna_packed());
                            produced_bases += len;
                        }
                    }

                    if (!produced.empty()) {
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            pending_bases += produced.bases;
                            batches.push_back(std::move(produced));
                            max_queue_depth = std::max(max_queue_depth, batches.size());
                            max_pending_bases = std::max(max_pending_bases, pending_bases);
                        }
                        ++worker_stats[worker_id].decompressed_chunks;
                        worker_stats[worker_id].decompressed_bytes += produced_bases;
                    }

                    if (!eof && !stop.load(std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        ready_streams.push_back(stream_id);
                    } else {
                        streams[stream_id].parser.reset();
                        std::string next_path;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (next_file < n_files) {
                                next_path = cfg.input_files[next_file++];
                                ++opened_streams;
                            } else {
                                --active_streams;
                            }
                        }
                        if (!next_path.empty()) {
                            GzInput inp(next_path);
                            auto next_parser = std::make_unique<AdaptivePackedFastqParser>(std::move(inp));
                            {
                                std::lock_guard<std::mutex> lk(q_mutex);
                                streams[stream_id].parser = std::move(next_parser);
                                ready_streams.push_back(stream_id);
                            }
                        }
                    }
                    q_cv.notify_all();
                    continue;
                }

                if (do_parse) {
                    ++worker_stats[worker_id].parse_tasks;
                    ++worker_stats[worker_id].parsed_chunks;
                    worker_stats[worker_id].records += batch.records.size();
                    worker_stats[worker_id].bases += batch.bases;
                    for (const auto& rec : batch.records) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        PackedReadView chunk{
                            batch.words.data() + rec.word_offset,
                            rec.word_count,
                            rec.tail,
                            rec.len
                        };
#ifdef TUNA_PACKED_NT_PHASE1
                        extract_superkmers_from_packed_nt<k, m>(
                            chunk, partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#else
                        decode_packed_nt_to_kache(chunk, kache_buf);
                        extract_superkmers_from_kache<k, m>(
                            kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                            local_kmers, local_superkmers, flush_fn);
#endif
                        ++local_seqs;
                    }
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            worker_stats[worker_id].kmers += local_kmers;
            worker_stats[worker_id].superkmers += local_superkmers;
        } catch (...) {
            record_error();
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);
    if (queue_stats) {
        phase1_print_adaptive_stats(
            "mem_adaptive_packed_gz_fastq", worker_stats, max_queue_depth,
            max_pending_bases, low_batch_bases, high_batch_bases, opened_streams);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Producer-consumer harness for a single gz file ──────────────────────────
//
// When there is only one input file and it is gzip-compressed, the standard
// work-stealing harness wastes all but one thread (min(n_threads,1)=1).
// This harness instead:
//   • 1 producer thread: gz decompression + split_actg → batches of ACTG chunks
//   • n_threads-1 consumer threads: extract_superkmers_from_actg per chunk
//
// The queue is bounded (MAX_QUEUE batches) for backpressure.  Consumers flush
// SuperkmerWriters to the shared bucket files under per-bucket mutexes.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_gz_pc(
    const Config&               cfg,
    const std::string&          gz_path,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      n_threads,          // ≥ 2 (1 producer + rest consumers)
    size_t                      write_budget_per_thread)
{
    using Batch = std::vector<std::string>;     // ACTG-only chunks pre-split by producer

    constexpr size_t MAX_QUEUE  = 32;           // max batches in flight
    constexpr size_t BATCH_SEQS = 512;          // sequences per batch

    const size_t n_parts     = cfg.num_partitions;
    const size_t n_consumers = n_threads - 1;

    std::deque<Batch>       queue;
    std::mutex              q_mutex;
    std::condition_variable q_cv;
    bool                    producer_done = false;
    std::exception_ptr      producer_error = nullptr;
    std::atomic<bool>       stop{false};
    std::exception_ptr      consumer_error = nullptr;
    std::mutex              consumer_error_mutex;

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

    // Producer: decompress gz via GzInput, use helicase SIMD parser to deliver
    // ACTG-only chunks, accumulate into batches and push to the queue.
    auto producer_fn = [&]() {
        auto feed = [&](auto& parser) {
            while (true) {
                Batch batch;
                size_t chunk_count = 0;
                while (!stop.load(std::memory_order_relaxed)
                       && chunk_count < BATCH_SEQS
                       && parser.next()) {
                    auto [ptr, len] = parser.get_dna_raw();
                    batch.emplace_back(ptr, len);
                    ++chunk_count;
                }
                if (stop.load(std::memory_order_relaxed)) break;
                if (batch.empty()) break;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{ return queue.size() < MAX_QUEUE; });
                    queue.push_back(std::move(batch));
                }
                q_cv.notify_one();
            }
            { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
            q_cv.notify_all();
        };
        try {
            GzInput inp(gz_path);
            if (inp.first_byte() == '@') {
                helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            } else {
                helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_error = std::current_exception();
                producer_done = true;
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    // Consumer: pull batches from the queue and extract superkmers.
    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                    if (queue.empty()) break;       // producer_done && queue empty → done
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();                  // wake producer if it was waiting for space

                for (const auto& chunk : batch) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    extract_superkmers_from_actg<k, m>(
                        chunk.data(), chunk.size(), partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(consumer_error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// Single gzipped FASTQ path with KMC-style raw chunking: the producer only
// decompresses and cuts at record boundaries; consumers parse/decode chunks.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_single_gz_fastq_raw_pc(
    const Config&               cfg,
    const std::string&          gz_path,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = RawFastqChunk;

    constexpr size_t MAX_QUEUE = 8;
#ifdef TUNA_FASTQ_CHUNK_MB
    constexpr size_t TARGET_CHUNK_BYTES = static_cast<size_t>(TUNA_FASTQ_CHUNK_MB) << 20;
#else
    constexpr size_t TARGET_CHUNK_BYTES = 4u << 20;
#endif

    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_consumers = std::max<size_t>(1, n_threads - 1);

    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    bool producer_done = false;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    size_t max_queue_depth = 0;
    Phase1QueueThreadStats producer_stats;
    std::vector<Phase1QueueThreadStats> consumer_stats(n_consumers);

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto producer_fn = [&]() {
        try {
            RawGzFastqChunker chunker(gz_path, TARGET_CHUNK_BYTES);
            Batch batch;
            while (!stop.load(std::memory_order_relaxed) && chunker.next(batch)) {
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) producer_stats.wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;
                    if (queue_stats) {
                        ++producer_stats.batches;
                        producer_stats.bases += batch.data.size();
                        max_queue_depth = std::max(max_queue_depth, queue.size() + 1);
                    }
                    queue.push_back(std::move(batch));
                }
                q_cv.notify_one();
                batch = Batch{};
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
        {
            std::lock_guard<std::mutex> lk(q_mutex);
            producer_done = true;
        }
        q_cv.notify_all();
    };

    auto consumer_fn = [&](size_t consumer_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return !queue.empty() || producer_done || stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) consumer_stats[consumer_id].wait_ns += phase1_now_ns() - wait_start;
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();

                if (batch.data.empty()) continue;
                if (queue_stats) {
                    ++consumer_stats[consumer_id].batches;
                    consumer_stats[consumer_id].bases += batch.data.size();
                }
                helicase::FastqParser<HELICASE_ACTG_PACKED, helicase::SliceInput> parser(
                    batch.data.data(), batch.data.size());
                while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                    const size_t len = parser.get_dna_len();
                    if (len < k) continue;
#ifdef TUNA_PACKED_NT_PHASE1
                    auto dna = parser.get_dna_packed();
                    auto [words, tail] = dna.bits();
                    PackedReadView rec{words.data(), words.size(), tail, dna.len()};
                    extract_superkmers_from_packed_nt<k, m>(
                        rec, partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#else
                    decode_packed_nt_to_kache(parser.get_dna_packed(), kache_buf);
                    extract_superkmers_from_kache<k, m>(
                        kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#endif
                    ++local_seqs;
                    if (queue_stats) {
                        ++consumer_stats[consumer_id].records;
                    }
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            if (queue_stats) {
                consumer_stats[consumer_id].kmers += local_kmers;
                consumer_stats[consumer_id].superkmers += local_superkmers;
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(1 + n_consumers);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn, t);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);
    if (queue_stats) {
        std::cerr << "[phase1-queue] single_gz_fastq_raw"
                  << " max_queue_depth=" << max_queue_depth
                  << " producer_batches=" << producer_stats.batches
                  << " producer_bytes=" << producer_stats.bases
                  << " producer_wait_s=" << phase1_ns_to_s(producer_stats.wait_ns)
                  << " consumers=" << n_consumers << "\n";
        phase1_print_thread_stats("single_gz_fastq_raw_consumer", consumer_stats);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// In-memory variant of the raw single-gzip FASTQ path.  The producer only
// decompresses and cuts at record boundaries; consumers parse/decode chunks and
// flush superkmers directly to per-partition memory buffers.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_single_gz_fastq_raw_pc(
    const Config&             cfg,
    const std::string&        gz_path,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    using Batch = RawFastqChunk;

    constexpr size_t MAX_QUEUE = 8;
#ifdef TUNA_FASTQ_CHUNK_MB
    constexpr size_t TARGET_CHUNK_BYTES = static_cast<size_t>(TUNA_FASTQ_CHUNK_MB) << 20;
#else
    constexpr size_t TARGET_CHUNK_BYTES = 4u << 20;
#endif

    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_consumers = std::max<size_t>(1, n_threads - 1);

    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    bool producer_done = false;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    size_t max_queue_depth = 0;
    Phase1QueueThreadStats producer_stats;
    std::vector<Phase1QueueThreadStats> consumer_stats(n_consumers);

    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto producer_fn = [&]() {
        try {
            RawGzFastqChunker chunker(gz_path, TARGET_CHUNK_BYTES);
            Batch batch;
            while (!stop.load(std::memory_order_relaxed) && chunker.next(batch)) {
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) producer_stats.wait_ns += phase1_now_ns() - wait_start;
                    if (stop.load(std::memory_order_relaxed)) break;
                    if (queue_stats) {
                        ++producer_stats.batches;
                        producer_stats.bases += batch.data.size();
                        max_queue_depth = std::max(max_queue_depth, queue.size() + 1);
                    }
                    queue.push_back(std::move(batch));
                }
                q_cv.notify_one();
                batch = Batch{};
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
        {
            std::lock_guard<std::mutex> lk(q_mutex);
            producer_done = true;
        }
        q_cv.notify_all();
    };

    auto consumer_fn = [&](size_t consumer_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return !queue.empty() || producer_done || stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) consumer_stats[consumer_id].wait_ns += phase1_now_ns() - wait_start;
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();

                if (batch.data.empty()) continue;
                if (queue_stats) {
                    ++consumer_stats[consumer_id].batches;
                    consumer_stats[consumer_id].bases += batch.data.size();
                }
                helicase::FastqParser<HELICASE_ACTG_PACKED, helicase::SliceInput> parser(
                    batch.data.data(), batch.data.size());
                while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                    const size_t len = parser.get_dna_len();
                    if (len < k) continue;
#ifdef TUNA_PACKED_NT_PHASE1
                    auto dna = parser.get_dna_packed();
                    auto [words, tail] = dna.bits();
                    PackedReadView rec{words.data(), words.size(), tail, dna.len()};
                    extract_superkmers_from_packed_nt<k, m>(
                        rec, partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#else
                    decode_packed_nt_to_kache(parser.get_dna_packed(), kache_buf);
                    extract_superkmers_from_kache<k, m>(
                        kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#endif
                    ++local_seqs;
                    if (queue_stats)
                        ++consumer_stats[consumer_id].records;
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            if (queue_stats) {
                consumer_stats[consumer_id].kmers += local_kmers;
                consumer_stats[consumer_id].superkmers += local_superkmers;
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(1 + n_consumers);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn, t);
    for (auto& th : threads) th.join();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);
    if (queue_stats) {
        std::cerr << "[phase1-queue] mem_single_gz_fastq_raw"
                  << " max_queue_depth=" << max_queue_depth
                  << " producer_batches=" << producer_stats.batches
                  << " producer_bytes=" << producer_stats.bases
                  << " producer_wait_s=" << phase1_ns_to_s(producer_stats.wait_ns)
                  << " consumers=" << n_consumers << "\n";
        phase1_print_thread_stats("mem_single_gz_fastq_raw_consumer", consumer_stats);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// Multiple gzipped inputs need finer-grained scheduling than file-level work
// stealing: each producer owns one gz stream at a time and consumers partition
// copied ACTG chunks from a shared queue.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_multi_gz_pc(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_files = cfg.input_files.size();
    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_producers = phase1_gz_producer_threads(n_files, n_threads);
    const size_t n_consumers = std::max<size_t>(1, n_threads - n_producers);

    std::atomic<size_t> next_file{0};
    std::atomic<size_t> active_producers{n_producers};
    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    size_t max_queue_depth = 0;
    std::vector<Phase1QueueThreadStats> producer_stats(n_producers);
    std::vector<Phase1QueueThreadStats> consumer_stats(n_consumers);

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases, size_t producer_id) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (queue_stats) producer_stats[producer_id].wait_ns += phase1_now_ns() - wait_start;
            if (stop.load(std::memory_order_relaxed)) return;
            if (queue_stats) {
                ++producer_stats[producer_id].batches;
                producer_stats[producer_id].records += batch.records.size();
                producer_stats[producer_id].bases += batch_bases;
                max_queue_depth = std::max(max_queue_depth, queue.size() + 1);
            }
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
        batch_bases = 0;
    };

    auto producer_done = [&]() {
        active_producers.fetch_sub(1, std::memory_order_relaxed);
        q_cv.notify_all();
    };

    auto producer_fn = [&](size_t producer_id) {
        try {
            while (!stop.load(std::memory_order_relaxed)) {
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;
                const auto& input_path = cfg.input_files[fi];

                Batch batch;
                batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                size_t batch_bases = 0;
                GzInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases, producer_id);
                        }
                    }
                } else {
                    helicase::FastaParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases, producer_id);
                        }
                    }
                }
                if (!stop.load(std::memory_order_relaxed))
                    push_batch(batch, batch_bases, producer_id);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
        producer_done();
    };

    auto consumer_fn = [&](size_t consumer_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return !queue.empty() ||
                               active_producers.load(std::memory_order_relaxed) == 0 ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) consumer_stats[consumer_id].wait_ns += phase1_now_ns() - wait_start;
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();
                if (queue_stats) {
                    ++consumer_stats[consumer_id].batches;
                    consumer_stats[consumer_id].records += batch.records.size();
                    for (const auto& rec : batch.records)
                        consumer_stats[consumer_id].bases += rec.len;
                }

                for (const auto& rec : batch.records) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    PackedReadView chunk{
                        batch.words.data() + rec.word_offset,
                        rec.word_count,
                        rec.tail,
                        rec.len
                    };
#ifdef TUNA_PACKED_NT_PHASE1
                    extract_superkmers_from_packed_nt<k, m>(
                        chunk, partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#else
                    decode_packed_nt_to_kache(chunk, kache_buf);
                    extract_superkmers_from_kache<k, m>(
                        kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#endif
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            if (queue_stats) {
                consumer_stats[consumer_id].kmers += local_kmers;
                consumer_stats[consumer_id].superkmers += local_superkmers;
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_producers + n_consumers);
    for (size_t t = 0; t < n_producers; ++t)
        threads.emplace_back(producer_fn, t);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn, t);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);
    if (queue_stats) {
        std::cerr << "[phase1-queue] multi_gz"
                  << " producers=" << n_producers
                  << " consumers=" << n_consumers
                  << " max_queue_depth=" << max_queue_depth << "\n";
        phase1_print_thread_stats("multi_gz_producer", producer_stats);
        phase1_print_thread_stats("multi_gz_consumer", consumer_stats);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// One or a few large plain inputs need finer-grained scheduling than file-level
// work stealing.  A producer mmaps each file, publishes string_views into that
// mapping, and waits for consumers to drain them before advancing to the next
// file so the views never outlive their backing mmap.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_plain_pc(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = std::vector<std::string_view>;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 8192;
    constexpr size_t TARGET_BATCH_BASES = 1u << 20;
    constexpr size_t TARGET_CHUNK_BASES = 1u << 20;

    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_consumers = n_threads - 1;

    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    bool producer_done = false;
    size_t active_batches = 0;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (stop.load(std::memory_order_relaxed)) return;
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch_bases = 0;
    };

    auto wait_until_consumed = [&]() {
        std::unique_lock<std::mutex> lk(q_mutex);
        q_cv.wait(lk, [&]{
            return (queue.empty() && active_batches == 0) ||
                   stop.load(std::memory_order_relaxed);
        });
    };

    auto producer_fn = [&]() {
        auto feed = [&](auto& parser) {
            Batch batch;
            size_t batch_bases = 0;
            auto add_chunk = [&](const char* ptr, size_t len) {
                if (len < k) return;
                const size_t step = TARGET_CHUNK_BASES > k
                    ? TARGET_CHUNK_BASES - (k - 1)
                    : TARGET_CHUNK_BASES;
                for (size_t off = 0; off < len; ) {
                    const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                    batch.emplace_back(ptr + off, sub_len);
                    batch_bases += sub_len;
                    if (batch_bases >= TARGET_BATCH_BASES ||
                        batch.size() >= MAX_BATCH_ITEMS)
                        push_batch(batch, batch_bases);
                    if (off + sub_len >= len) break;
                    off += step;
                }
            };

            while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                auto [ptr, len] = parser.get_dna_raw();
                add_chunk(ptr, len);
            }
            if (!stop.load(std::memory_order_relaxed))
                push_batch(batch, batch_bases);
            wait_until_consumed();
        };

        try {
            for (const auto& input_path : cfg.input_files) {
                if (stop.load(std::memory_order_relaxed)) break;
                helicase::MmapInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    feed(p);
                } else {
                    helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    feed(p);
                }
            }
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_done = true;
            }
            q_cv.notify_all();
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_error = std::current_exception();
                producer_done = true;
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return !queue.empty() || producer_done ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                    ++active_batches;
                }
                q_cv.notify_all();

                for (const auto chunk : batch) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    extract_superkmers_from_actg<k, m>(
                        chunk.data(), chunk.size(), partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                }

                {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    --active_batches;
                }
                q_cv.notify_all();
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// File-level work stealing: each worker atomically claims the next file and runs
// helicase + superkmer extraction on its own thread.
// Load balancing is coarse (file granularity) but negligible for equal-size files.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);

    if (n_files >= 1 && n_threads_req > 1) {
        bool all_gz = true;
        bool all_plain = true;
        for (const auto& f : cfg.input_files) {
            if (!(f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)) {
                all_gz = false;
            } else {
                all_plain = false;
            }
        }
        if (cfg.phase1_adaptive && all_gz) {
            if (phase1_adaptive_all_gz_fastq(cfg, all_gz)) {
                if (n_files == 1)
                    return partition_kmers_single_gz_fastq_raw_pc<k, m>(
                        cfg, cfg.input_files[0], buckets, partition_fn, write_budget_per_thread);
                return partition_kmers_adaptive_packed_gz_fastq_pc<k, m>(
                    cfg, buckets, partition_fn, write_budget_per_thread);
            }
            if (!cfg.hide_progress)
                std::cerr << "tuna: note: -p1-adaptive currently supports only gz FASTQ; using default phase 1\n";
        }
        if (all_gz && n_files == 1 && gz_first_byte(cfg.input_files[0]) == '@')
            return partition_kmers_single_gz_fastq_raw_pc<k, m>(
                cfg, cfg.input_files[0], buckets, partition_fn, write_budget_per_thread);
        if (all_gz)
            return partition_kmers_multi_gz_pc<k, m>(
                cfg, buckets, partition_fn, write_budget_per_thread);
        if (all_plain && n_files <= n_threads_req)
            return partition_kmers_plain_pc<k, m>(
                cfg, buckets, partition_fn, write_budget_per_thread);
    }

    const size_t n_threads = std::min(n_threads_req, n_files);
    const size_t n_parts   = cfg.num_partitions;

    std::atomic<size_t>   next_file{0};
    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads_req / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread),
        cfg.lz4_buckets);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            SeqSource            source;
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (worker_error) std::rethrow_exception(worker_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Parallel harness (in-memory sinks) ───────────────────────────────────────
//
// Same logic as partition_kmers_impl but writes packed superkmers into
// per-partition std::string buffers instead of disk files.  Used by the
// streaming pipeline to avoid the Phase 1 disk write + Phase 2 mmap round-trip.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_multi_gz_pc(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_files = cfg.input_files.size();
    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_producers = phase1_gz_producer_threads(n_files, n_threads);
    const size_t n_consumers = std::max<size_t>(1, n_threads - n_producers);

    std::atomic<size_t> next_file{0};
    std::atomic<size_t> active_producers{n_producers};
    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;
    const bool queue_stats = phase1_queue_stats_enabled();
    size_t max_queue_depth = 0;
    std::vector<Phase1QueueThreadStats> producer_stats(n_producers);
    std::vector<Phase1QueueThreadStats> consumer_stats(n_consumers);

    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases, size_t producer_id) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (queue_stats) producer_stats[producer_id].wait_ns += phase1_now_ns() - wait_start;
            if (stop.load(std::memory_order_relaxed)) return;
            if (queue_stats) {
                ++producer_stats[producer_id].batches;
                producer_stats[producer_id].records += batch.records.size();
                producer_stats[producer_id].bases += batch_bases;
                max_queue_depth = std::max(max_queue_depth, queue.size() + 1);
            }
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
        batch_bases = 0;
    };

    auto producer_done = [&]() {
        active_producers.fetch_sub(1, std::memory_order_relaxed);
        q_cv.notify_all();
    };

    auto producer_fn = [&](size_t producer_id) {
        try {
            while (!stop.load(std::memory_order_relaxed)) {
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;
                const auto& input_path = cfg.input_files[fi];

                Batch batch;
                batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                size_t batch_bases = 0;
                GzInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases, producer_id);
                        }
                    }
                } else {
                    helicase::FastaParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases, producer_id);
                        }
                    }
                }
                if (!stop.load(std::memory_order_relaxed))
                    push_batch(batch, batch_bases, producer_id);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
        producer_done();
    };

    auto consumer_fn = [&](size_t consumer_id) {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    const uint64_t wait_start = queue_stats ? phase1_now_ns() : 0;
                    q_cv.wait(lk, [&]{
                        return !queue.empty() ||
                               active_producers.load(std::memory_order_relaxed) == 0 ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue_stats) consumer_stats[consumer_id].wait_ns += phase1_now_ns() - wait_start;
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();
                if (queue_stats) {
                    ++consumer_stats[consumer_id].batches;
                    consumer_stats[consumer_id].records += batch.records.size();
                    for (const auto& rec : batch.records)
                        consumer_stats[consumer_id].bases += rec.len;
                }

                for (const auto& rec : batch.records) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    PackedReadView chunk{
                        batch.words.data() + rec.word_offset,
                        rec.word_count,
                        rec.tail,
                        rec.len
                    };
#ifdef TUNA_PACKED_NT_PHASE1
                    extract_superkmers_from_packed_nt<k, m>(
                        chunk, partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#else
                    decode_packed_nt_to_kache(chunk, kache_buf);
                    extract_superkmers_from_kache<k, m>(
                        kache_buf.data(), kache_buf.size(), partition_fn, min_it, writers,
                        local_kmers, local_superkmers, flush_fn);
#endif
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
            if (queue_stats) {
                consumer_stats[consumer_id].kmers += local_kmers;
                consumer_stats[consumer_id].superkmers += local_superkmers;
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_producers + n_consumers);
    for (size_t t = 0; t < n_producers; ++t)
        threads.emplace_back(producer_fn, t);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn, t);
    for (auto& th : threads) th.join();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);
    if (queue_stats) {
        std::cerr << "[phase1-queue] mem_multi_gz"
                  << " producers=" << n_producers
                  << " consumers=" << n_consumers
                  << " max_queue_depth=" << max_queue_depth << "\n";
        phase1_print_thread_stats("mem_multi_gz_producer", producer_stats);
        phase1_print_thread_stats("mem_multi_gz_consumer", consumer_stats);
    }

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_impl(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts      = cfg.num_partitions;

    // Single file + multiple threads → producer-consumer (reuse the same
    // consumers but flush to bufs instead of ofstreams).
    const size_t n_threads = (n_files == 1 && n_threads_req > 1)
        ? n_threads_req   // single-file producer-consumer: all threads participate
        : std::min(n_threads_req, n_files);

    if (n_files >= 1 && n_threads_req > 1) {
        bool all_gz = true;
        for (const auto& input_path : cfg.input_files) {
            if (!(input_path.size() > 3 &&
                  input_path.compare(input_path.size() - 3, 3, ".gz") == 0)) {
                all_gz = false;
                break;
            }
        }
        if (cfg.phase1_adaptive && all_gz) {
            if (phase1_adaptive_all_gz_fastq(cfg, all_gz)) {
                if (n_files == 1)
                    return partition_kmers_mem_single_gz_fastq_raw_pc<k, m>(
                        cfg, cfg.input_files[0], bufs, partition_fn);
                return partition_kmers_mem_adaptive_packed_gz_fastq_pc<k, m>(
                    cfg, bufs, partition_fn);
            }
            if (!cfg.hide_progress)
                std::cerr << "tuna: note: -p1-adaptive currently supports only gz FASTQ; using default phase 1\n";
        }
        if (all_gz && n_files == 1 && gz_first_byte(cfg.input_files[0]) == '@')
            return partition_kmers_mem_single_gz_fastq_raw_pc<k, m>(
                cfg, cfg.input_files[0], bufs, partition_fn);
        if (all_gz)
            return partition_kmers_mem_multi_gz_pc<k, m>(cfg, bufs, partition_fn);
    }

    // Single file: producer-consumer variant using in-memory sinks.
    if (n_files == 1 && n_threads_req > 1) {
        const auto& input_path = cfg.input_files[0];
        const bool is_gz = input_path.size() > 3 &&
                           input_path.compare(input_path.size() - 3, 3, ".gz") == 0;

        if (!is_gz) {
            using Batch = std::vector<std::string_view>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            size_t                  active_batches = 0;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };

                auto wait_until_consumed = [&]() {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return (queue.empty() && active_batches == 0) ||
                               stop.load(std::memory_order_relaxed);
                    });
                };

                auto feed = [&](auto& parser) {
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };

                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        auto [ptr, len] = parser.get_dna_raw();
                        add_chunk(ptr, len);
                    }
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_done = true;
                    }
                    q_cv.notify_all();
                    wait_until_consumed();
                };

                try {
                    helicase::MmapInput inp(input_path);
                    if (inp.first_byte() == '@') {
                        helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                        feed(p);
                    } else {
                        helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                        feed(p);
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                            if (queue.empty()) break;
                            batch = std::move(queue.front());
                            queue.pop_front();
                            ++active_batches;
                        }
                        q_cv.notify_all();
                        for (const auto chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            --active_batches;
                        }
                        q_cv.notify_all();
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
        }

            using Batch = std::vector<std::string>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };
                try {
                    SeqSource source;
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };
                    source.process(input_path, [&](const char* ptr, size_t len) {
                        if (stop.load(std::memory_order_relaxed)) return;
                        add_chunk(ptr, len);
                    });
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
                    q_cv.notify_all();
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        { std::unique_lock<std::mutex> lk(q_mutex);
                          q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                          if (queue.empty()) break;
                          batch = std::move(queue.front());
                          queue.pop_front(); }
                        q_cv.notify_one();
                        for (const auto& chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
    }

    // Few large plain files: split parsed records/chunks across all worker
    // threads instead of assigning whole files to workers.  The producer keeps
    // each mmap-backed parser alive until consumers have drained all views for
    // that file.
    if (n_files > 1 && n_files <= n_threads_req && n_threads_req > 1) {
        bool all_plain = true;
        for (const auto& input_path : cfg.input_files) {
            if (input_path.size() > 3 &&
                input_path.compare(input_path.size() - 3, 3, ".gz") == 0) {
                all_plain = false;
                break;
            }
        }

        if (all_plain) {
            using Batch = std::vector<std::string_view>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads_req - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            size_t                  active_batches = 0;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };

                auto wait_until_consumed = [&]() {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return (queue.empty() && active_batches == 0) ||
                               stop.load(std::memory_order_relaxed);
                    });
                };

                auto feed = [&](auto& parser) {
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };

                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        auto [ptr, len] = parser.get_dna_raw();
                        add_chunk(ptr, len);
                    }
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    wait_until_consumed();
                };

                try {
                    for (const auto& input_path : cfg.input_files) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        helicase::MmapInput inp(input_path);
                        if (inp.first_byte() == '@') {
                            helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                            feed(p);
                        } else {
                            helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                            feed(p);
                        }
                    }
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_done = true;
                    }
                    q_cv.notify_all();
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                            if (queue.empty()) break;
                            batch = std::move(queue.front());
                            queue.pop_front();
                            ++active_batches;
                        }
                        q_cv.notify_all();
                        for (const auto chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            --active_batches;
                        }
                        q_cv.notify_all();
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads_req);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
        }
    }

    // Multi-file: file-level work-stealing.
    std::atomic<size_t>     next_file{0};
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            SeqSource            source;
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);
    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Public ───────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers(
    const Config&               cfg,
    std::vector<SuperkmerBucketFile>& buckets,
    size_t                      write_budget_per_thread)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    return partition_kmers_impl<k, m>(cfg, buckets,
        [mask](uint64_t h) -> size_t { return h & mask; },
        write_budget_per_thread);
}

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_mem(
    const Config&             cfg,
    std::vector<std::string>& bufs)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    return partition_kmers_mem_impl<k, m>(cfg, bufs,
        [mask](uint64_t h) -> size_t { return h & mask; });
}
