
#ifndef KACHE_HASH_UTILITY_HPP
#define KACHE_HASH_UTILITY_HPP



#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <cassert>


namespace kache_hash
{


// Returns pointer to a memory-allocation for `size` elements of type `T_`.
template <typename T_>
inline T_* allocate(const std::size_t size)
{
    return static_cast<T_*>(std::malloc(size * sizeof(T_)));
}


// Returns pointer to a memory-allocation for `size` elements of type `T_`,
// whose alignment is specified is `alignment`.
template <typename T_>
inline T_* aligned_allocate(const std::size_t size, const std::size_t alignment = 8)
{
    const auto bytes = ((size * sizeof(T_) + alignment - 1) / alignment) * alignment;
    return static_cast<T_*>(std::aligned_alloc(alignment, bytes));
}


// Returns pointer to a memory-reallocation for `size` elements of type `T_`.
template <typename T_>
inline T_* reallocate(T_* const ptr, const std::size_t size)
{
    return static_cast<T_*>(std::realloc(ptr, size * sizeof(T_)));
}


// Deallocates the pointer `ptr`, allocated with `allocate`.
template <typename T_>
inline void deallocate(T_* const ptr)
{
    std::free(ptr);
}


// Serializes `x` to the stream `os`.
template <typename T_>
inline void serialize(const T_ x, std::ofstream& os)
{
    os.write(reinterpret_cast<const char*>(&x), sizeof(x));
}


// Serializes the array `A` of size `n` to the stream `os`.
template <typename T_>
inline void serialize(const T_* const A, const std::size_t n, std::ofstream& os)
{
    os.write(reinterpret_cast<const char*>(A), n * sizeof(T_));
}


// Deserializes `x` from the stream `is`.
template <typename T_>
inline void deserialize(T_& x, std::ifstream& is)
{
    is.read(reinterpret_cast<char*>(&x), sizeof(x));
}


// Deserializes the array `A` of size `n` from the stream `is`.
template <typename T_>
inline void deserialize(T_* const A, std::size_t n, std::ifstream& is)
{
    is.read(reinterpret_cast<char*>(A), n * sizeof(T_));
}


// Returns the smallest power of 2 at least as large as `x`. `x` must be in
// `[1, 2^63]`.
inline constexpr std::size_t ceil_pow_2(std::size_t x)
{
    assert(x > 0 && x <= (1lu << 63));

    if((x & (x - 1)) == 0)
        return x;

    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);

    return x + 1;
}


// Returns the bits in the indices `[l, h]` (0-indexed) of `x`.
inline constexpr uint64_t bits(uint64_t x, uint64_t l, uint64_t h)
{
    assert(l <= h && h < 64);
    return (x << (63 - h)) >> ((63 - h) + l);
}


// Fields for profiling time.
typedef std::chrono::high_resolution_clock::time_point time_point_t;
constexpr static auto now = std::chrono::high_resolution_clock::now;
constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


// Wrapper class for a data element of type `T_` to ensure that in a linear
// collection of `T_`'s, each element is aligned to a cache-line boundary.
template <typename T_>
class alignas(L1_CACHE_LINE_SIZE)
    Padded
{
private:

    T_ data_;


public:

    Padded()
    {}

    Padded(const T_& data):
      data_(data)
    {}

    Padded(T_&& data):
        data_(std::move(data))
    {}

    T_& unwrap() { return data_; }

    const T_& unwrap() const { return data_; }
};


}



#endif
