#pragma once

// Thread-safe KFF output writer.
//
// Wraps a Kff_file + Section_Raw; multiple threads can call write_batch
// concurrently — serialised by an internal mutex.
//
// Encoding: A=0, C=1, G=2, T=3 (standard 2-bit, MSB-first packing).
// Flags: canonical=true (tuna reports the lex-min k-mer), unique=true.
// data_size=4: each k-mer carries a uint32_t count stored big-endian.

#include "kff_io.hpp"

#include <cstdint>
#include <mutex>
#include <string>

class KffOutput {
public:
    explicit KffOutput(const std::string& path, uint16_t k)
        : file_(path, "w"), k_(k), kbytes_((k + 3) / 4)
    {
        file_.write_encoding(0, 1, 2, 3);  // A=0 C=1 G=2 T=3
        file_.set_uniqueness(true);
        file_.set_canonicity(true);
        file_.write_metadata(0, nullptr);

        Section_GV gv(&file_);
        gv.write_var("k", k_);
        gv.write_var("max", 1);        // one k-mer per block
        gv.write_var("data_size", 4);  // uint32_t count per k-mer
        gv.close();

        raw_ = new Section_Raw(&file_);
    }

    ~KffOutput() { close(); }

    // Write n_kmers packed k-mers.
    // packed: flat byte array, kbytes_ bytes per k-mer (2-bit, MSB-first).
    // counts: one uint32_t per k-mer (written big-endian).
    void write_batch(const uint8_t* packed, const uint32_t* counts, size_t n_kmers) {
        uint8_t cbuf[4];
        std::lock_guard<std::mutex> g(mutex_);
        for (size_t i = 0; i < n_kmers; ++i) {
            const uint32_t c = counts[i];
            cbuf[0] = (c >> 24) & 0xFF;
            cbuf[1] = (c >> 16) & 0xFF;
            cbuf[2] = (c >>  8) & 0xFF;
            cbuf[3] =  c        & 0xFF;
            raw_->write_compacted_sequence(
                // kff API takes non-const pointer but only reads the data
                const_cast<uint8_t*>(packed + i * kbytes_),
                static_cast<uint64_t>(k_),
                cbuf);
        }
    }

    void close() {
        if (raw_) {
            raw_->close();
            delete raw_;
            raw_ = nullptr;
        }
    }

private:
    Kff_file     file_;
    Section_Raw* raw_ = nullptr;
    std::mutex   mutex_;
    uint16_t     k_;
    size_t       kbytes_;
};
