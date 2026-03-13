
#ifndef KACHE_HASH_CONCURRENT_HASH_TABLE_HPP
#define KACHE_HASH_CONCURRENT_HASH_TABLE_HPP



#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "xxHash/xxhash.h"

#include <cstddef>
#include <vector>
#include <cstring>
#include <optional>
#include <cstdlib>
#include <algorithm>
#include <cmath>


namespace kache_hash
{


// Hasher for objects of type `T_`.
template <typename T_>
struct Hash
{
    uint64_t operator()(const T_& key) const { return XXH3_64bits(&key, sizeof(key)); }
};


// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table
{
    static_assert(sizeof(T_key_) <= 8, "Key for concurrent hash table must have size at most 8 bytes for atomic-read guarantees to hold.");

    class Iterator;
    friend class Iterator;

private:

    struct Key_Val_Pair
    {
        T_key_ key;
        T_val_ val;
    };

    static constexpr double lf_default = 0.8;   // Default maximum load-factor supported.

    T_key_ empty_key_;  // The empty key; currently it's set to all 1-bits.

    const T_hasher_ hash;   // The hasher object.

    const std::size_t capacity_;    // True capacity of the table; adjusted to be a power of 2.
    const std::size_t idx_wrapper_mask; // Bitmask to wrap indexing into the table.

    Key_Val_Pair* const T;  // The flat table of key-value collection.

    mutable std::vector<Spin_Lock> lock;    // Locks for atomic storing of key-value entries.
    // TODO: experiment moving the spin-lock to the key-value structure for seemingly better cache-behavior.

    // Maps the hash value `h` to an index into the table.
    std::size_t hash_to_idx(std::size_t h) const { return h & idx_wrapper_mask; }

    // Returns the next (wrapped) index for `i`.
    std::size_t next_index(std::size_t i) const { return hash_to_idx(i + 1); }


public:

    // Constructs a concurrent hash table to support upto `max_n` elements, with
    // a maximum load-factor of `load_factor`. The object `hasher` is used to
    // hash the keys in the table.
    Concurrent_Hash_Table(std::size_t max_n, double load_factor = lf_default, T_hasher_ hasher = T_hasher_());

    ~Concurrent_Hash_Table() { deallocate(T); }

    // Returns the capacity of the hash table.
    std::size_t capacity() const { return capacity_; }

    // Clears the hash table.
    void clear();

    // Returns the size of the table. It costs linear in the size of the table.
    std::size_t size() const;

    // Inserts the key `key` with value `val` into the table. Returns `nullptr`
    // iff the insertion succeeds, i.e. the key was absent in the table. Returns
    // the address of the associated value otherwise.
    const T_val_* insert(T_key_ key, T_val_ val);

    // Inserts the key `key` with value `val` into the table if it is absent and
    // returns `null_opt`. Otherwise updates the associated value to `val` and
    // returns the old value.
    std::optional<T_val_> upsert(T_key_ key, T_val_ val);

    // Searches for `key` in the table and returns the address of the value
    // associated to it iff it is found. Returns `nullptr` otherwise.
    const T_val_* find(T_key_ key) const;

    // Finds `key` and applies `f` to its value in-place, returning the old
    // value. Returns `null_opt` iff the key is absent.
    template <typename F>
    std::optional<T_val_> find_and_update(T_key_ key, const F& f);

    // Returns an iterator for the key-value pairs in the table.
    Iterator iterator() { return Iterator(*this, 1, 0); }

    // Returns an iterator for the key-value pairs in the table that belongs to
    // a group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator iterator(const std::size_t it_count, const std::size_t it_id) { return Iterator(*this, it_count, it_id); }
};


// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator
{
    friend class Concurrent_Hash_Table;

private:

    Concurrent_Hash_Table& M;   // The hash table to iterate on.
    std::size_t idx;    // Current slot-index the iterator is in.
    std::size_t end;    // End-index of the slot-range for the iterator.


    // Constructs an iterator for the hash table `M`, that would belong to a
    // group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id);


public:

    // Moves the iterator to the next key in the table. Iff some key is found
    // within its remaining range, the key and the associated value are put at
    // `key` and `val` respectively, and returns true.
    bool next(T_key_& key, T_val_& val);
};


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Concurrent_Hash_Table(const std::size_t max_n, const double load_factor, const T_hasher_ hasher):
      empty_key_()
    , hash(hasher)
    , capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / load_factor))))
    , idx_wrapper_mask(capacity_ - 1)
    , T(allocate<Key_Val_Pair>(capacity_))
    , lock(capacity_)
{
    std::memset(reinterpret_cast<char*>(&empty_key_), -1, sizeof(empty_key_));

    clear();
}

template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::clear()
{
    std::memset(static_cast<void*>(T), -1, capacity_ * sizeof(Key_Val_Pair));
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline std::size_t Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::size() const
{
    std::size_t sz = 0;
    for(std::size_t i = 0; i < capacity_; ++i)
        sz += (T[i].key != empty_key_);

    return sz;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val)
{
    bool success = false;

    for(std::size_t tried = 0, i = hash_to_idx(hash(key)); tried < capacity_; ++tried, i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return nullptr;
        }

        if(T[i].key == key)
            return &T[i].val;
    }

    return &T[0].val;   // Table full — silently treat as already present.
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline std::optional<T_val_> Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::upsert(const T_key_ key, const T_val_ val)
{
    bool success = false;

    for(std::size_t tried = 0, i = hash_to_idx(hash(key)); tried < capacity_; ++tried, i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return std::nullopt;
        }

        if(T[i].key == key)
        {
            lock[i].lock();
            const auto old_val = T[i].val;
            T[i].val = val;
            lock[i].unlock();

            return old_val;
        }
    }

    return T[0].val;    // Table full — placeholder.
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key) const
{
#ifndef NDEBUG
    std::size_t tried_slots = 0;
#endif

    const T_val_* val_ptr = nullptr;
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        if(T[i].key == key)
        {
            lock[i].lock(); // To ensure that some other thread is not updating this val atm—a stable value is guaranteed.
            val_ptr = &T[i].val;
            lock[i].unlock();
            break;
        }
        else if(T[i].key == empty_key_)
            break;
        else
            assert(++tried_slots <= capacity_);
    }

    return val_ptr;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <typename F>
inline std::optional<T_val_> Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find_and_update(const T_key_ key, const F& f)
{
    for(std::size_t tried = 0, i = hash_to_idx(hash(key)); tried < capacity_; ++tried, i = next_index(i))
    {
        if(T[i].key == key)
        {
            lock[i].lock();
            const auto old_val = T[i].val;
            T[i].val = f(old_val);
            lock[i].unlock();
            return old_val;
        }
        else if(T[i].key == empty_key_)
            return std::nullopt;
    }
    return std::nullopt;  // Table full — key not found.
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id):
      M(M)
{
    const auto range_sz = (M.capacity_ + it_count - 1) / it_count;
    idx = it_id * range_sz;
    end = std::min((it_id + 1) * range_sz, M.capacity_);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::next(T_key_& key, T_val_& val)
{
    for(; idx < end; idx++)
        if(M.T[idx].key != M.empty_key_)
        {
            key = M.T[idx].key, val = M.T[idx].val;
            idx++;
            return true;
        }

    return false;
}


}



#endif