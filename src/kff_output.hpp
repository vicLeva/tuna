#pragma once

// Thread-safe batched KFF output writer.
//
// Encoding: A=0, C=1, G=2, T=3. KFF left-pads a partial first byte.
// Flags: unique=true and canonical reflects Tuna's -b setting.
// Counts use the smallest lossless 1-4 byte big-endian width per raw section.
//
// With max=1, each raw-section block has a fixed representation:
//
//     ceil(k / 4) sequence bytes + data_size count bytes
//
// There is no per-block k-mer-count field. Build complete interleaved batches
// in the counting threads and pass them directly to Kff_file::write(), avoiding
// two API calls and two buffered copies for every individual k-mer.

#include "kff_io.hpp"

#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>

class KffOutput {
public:
    explicit KffOutput(const std::string& path, uint16_t k, bool canonical)
        : file_(path, "w"), k_(k), kbytes_((k + 3) / 4)
    {
        // Tuna already submits megabyte-scale batches. Keep only the library's
        // small header buffer so record payloads go directly to the file
        // stream instead of being copied through a second 1 MB buffer.
        file_.max_buffer_size = file_.buffer_size;
        file_.write_encoding(0, 1, 2, 3);  // A=0 C=1 G=2 T=3
        file_.set_uniqueness(true);
        file_.set_canonicity(canonical);
        file_.write_metadata(0, nullptr);
    }

    ~KffOutput() { close(); }

    // records contains n_kmers complete raw-section records.
    void write_batch(
        const uint8_t* records,
        size_t n_kmers,
        uint8_t count_bytes)
    {
        if (n_kmers == 0) return;
        std::lock_guard<std::mutex> g(mutex_);
        if (count_bytes != count_bytes_) start_section(count_bytes);
        file_.write(records, n_kmers * (kbytes_ + count_bytes));
        raw_->nb_blocks += n_kmers;
    }

    void close() {
        if (closed_) return;
        if (!raw_) start_section(4);
        if (raw_) {
            raw_->close();
            delete raw_;
            raw_ = nullptr;
        }
        closed_ = true;
    }

private:
    void start_section(uint8_t count_bytes) {
        if (count_bytes < 1 || count_bytes > 4)
            throw std::invalid_argument("KFF count width must be in [1,4]");
        if (raw_) {
            raw_->close();
            delete raw_;
            raw_ = nullptr;
        }

        Section_GV gv(&file_);
        gv.write_var("k", k_);
        gv.write_var("max", 1);
        gv.write_var("data_size", count_bytes);
        gv.close();

        raw_ = new Section_Raw(&file_);
        count_bytes_ = count_bytes;
    }

    Kff_file     file_;
    Section_Raw* raw_ = nullptr;
    std::mutex   mutex_;
    uint16_t     k_;
    size_t       kbytes_;
    uint8_t      count_bytes_ = 0;
    bool         closed_ = false;
};
