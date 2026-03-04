
#ifndef KACHE_HASH_RW_LOCK_HPP
#define KACHE_HASH_RW_LOCK_HPP



#include "utility.hpp"

#include <cstddef>
#include <atomic>
#include <algorithm>


namespace kache_hash
{


// A readers-writer lock allowing concurrent read-only accesses and exclusive
// write access. It can support up-to `N` distinct workers.
template <std::size_t N = 128>
class RW_Lock
{

private:

    Padded<std::atomic_bool> reader[N]; // `reader[w]` denotes whether the `w`'th worker is making a read-access or not.
    std::atomic_flag writer;    // Whether a write access in ongoing or not.

public:

    constexpr RW_Lock();

    // Acquires exclusive ownership of the lock, blocking until no other worker
    // owns it in any form.
    void lock();

    // Releases exclusive ownership of the lock.
    void unlock();

    // Acquires shared ownership of the lock for worker `w_id`, blocking until
    // no worker owns it exclusively.
    void lock_shared(std::size_t w_id);

    // Releases shared ownership of the lock from worker `w_id`.
    void unlock_shared(std::size_t w_id);
};


template <std::size_t N>
inline constexpr RW_Lock<N>::RW_Lock():
      writer(false)
{
    std::for_each(reader, reader + N, [&](auto& r){ r.unwrap() = false; });
}


template <std::size_t N>
inline void RW_Lock<N>::lock()
{
    while(writer.test_and_set(std::memory_order_acquire))   // Wait until the ongoing write finishes, if any.
        writer.wait(true, std::memory_order_acquire);

    bool read;
    do  // Wait until ongoing reads finish.
    {
        read = false;
        std::for_each(reader, reader + N, [&](const auto& r){ read = read || r.unwrap(); });
    }
    while(read);
}


template <std::size_t N>
inline void RW_Lock<N>::unlock()
{
    assert(writer.test());
    writer.clear(std::memory_order_release);
    writer.notify_all();
}


template <std::size_t N>
inline void RW_Lock<N>::lock_shared(const std::size_t w_id)
{
    assert(w_id < N);
    auto& r = reader[w_id].unwrap();

    while(!r)   // Keep checking as long as it's not safe to commit to reading. This also makes the read-lock reentrant.
    {
        while(writer.test(std::memory_order_acq_rel))   // Some write is either (1) ongoing, in that case wait until it finishes; or
                                                        // (2) has signalled that it trying to start, in that case prioritize the write.
            writer.wait(true, std::memory_order_acq_rel);
        // Write(s) have finished.

        r = true;   // Declare intention to read. It is unsafe to exit right away, as some other writer may have started in between.
        if(writer.test(std::memory_order_acq_rel)) // A new writer has, observably, sneaked in in the mean time of read-declaration; unsafe to grab a read-lock.
            r = false;  // Give up the read-intention until the new writer finishes; unfair—continual stream of writers will starve readers.
                        // A strategy to prioritize reads can starve the writes, but writes typically being much more infrequent, read-starvation is preferable.
                        // A fair policy requires shared counters, incurring performance-hits.
    }

    // Committed to reading.
}


template <std::size_t N>
inline void RW_Lock<N>::unlock_shared(const std::size_t w_id)
{
    assert(w_id < N);   assert(reader[w_id].unwrap());
    reader[w_id].unwrap() = false;
}


}



#endif
