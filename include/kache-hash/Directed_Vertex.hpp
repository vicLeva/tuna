
#ifndef KACHE_HASH_DIRECTED_VERTEX_HPP
#define KACHE_HASH_DIRECTED_VERTEX_HPP



#include "DNA.hpp"
#include "Kmer.hpp"

#include <cstdint>


namespace kache_hash
{


// A class denoting an instance of a vertex. It's "directed" in the sense that the k-mer
// observed for the vertex is in a particular orientation — although a vertex `v` has an
// unambiguous canonical k-mer `v_hat`, the vertex can be observed in two different k-mer
// forms: `v_hat` and `{v_hat}_bar` — the class keeps track of the particular k-mer form
// observed for the vertex instance.
template <uint16_t k>
class Directed_Vertex
{
private:

    Kmer<k> kmer_;  // The observed k-mer for the vertex.
    Kmer<k> kmer_bar_;  // Reverse complement of the k-mer observed for the vertex.
    const Kmer<k>* kmer_hat_ptr;    // Pointer to the canonical form of the k-mer associated to the vertex.

    // Initializes the data of the class once the observed k-mer `kmer_` is set.
    void init();


public:

    // Constructs an empty vertex.
    Directed_Vertex()
    {}

    // Constructs a vertex observed for the k-mer `kmer`.
    Directed_Vertex(const Kmer<k>& kmer);

    // Copy constructs the vertex from `rhs`.
    Directed_Vertex(const Directed_Vertex<k>& rhs);

    // Assigns the vertex `rhs` to this one, and returns a constant reference to this object.
    const Directed_Vertex<k>& operator=(const Directed_Vertex<k>& rhs);

    // Returns `true` iff the k-mer observed for the vertex is in its canonical form.
    bool in_canonical_form() const;

    // Configures the vertex with the first k-mer from a super k-mer's binary
    // representation `super_kmer` that has `word_count` words. The super k-mer
    // is assumed to be MSB-aligned.
    void from_super_kmer(const uint64_t* super_kmer, std::size_t word_count);

    // Returns the observed k-mer for the vertex.
    const Kmer<k>& kmer() const;

    // Returns the reverse complement of the observed k-mer for the vertex.
    const Kmer<k>& kmer_bar() const;

    // Returns the canonical form of the vertex.
    const Kmer<k>& canonical() const;

    // Transforms this vertex to another by chopping off the first base from the associated
    // observed k-mer, and appending the nucleobase `b` to the end, i.e. effectively
    // rolling the associated k-mer by one base "forward".
    void roll_forward(DNA::Base b);

    // Returns a vertex formed by chopping off the first base from the observed
    // k-mer of this vertex, and appending the nucleobase `b` to the end, i.e.
    // effectively rolling the associated k-mer by one base "forward".
    const Directed_Vertex<k> roll_forward(DNA::Base b) const;

    // Transforms this vertex to another by chopping off the last base from the
    // associated observed k-mer, and appending the nucleobase `b` to the
    // beginning, i.e. effectively rolling the associated k-mer by one base
    // "backward".
    void roll_backward(DNA::Base b);

    // Returns a vertex formed by chopping off the last base from the observed
    // k-mer of this vertex, and appending the nucleobase `b` to the beginning,
    // i.e. effectively rolling the associated k-mer by one base "backward".
    const Directed_Vertex roll_backward(DNA::Base b) const;

    // Returns `true` iff this vertex and the vertex `v` are the same vertex, without the
    // directionality.
    bool is_same_vertex(const Directed_Vertex<k>& v) const;
};


template <uint16_t k>
inline void Directed_Vertex<k>::init()
{
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
}


template <uint16_t k>
inline Directed_Vertex<k>::Directed_Vertex(const Kmer<k>& kmer):
      kmer_(kmer)
{
    init();
}


template <uint16_t k>
inline Directed_Vertex<k>::Directed_Vertex(const Directed_Vertex<k>& rhs):
      kmer_(rhs.kmer_)
    , kmer_bar_(rhs.kmer_bar_)
    , kmer_hat_ptr(rhs.kmer_hat_ptr == &rhs.kmer_ ? &kmer_ : &kmer_bar_)
{}


template <uint16_t k>
inline const Directed_Vertex<k>& Directed_Vertex<k>::operator=(const Directed_Vertex<k>& rhs)
{
    kmer_ = rhs.kmer_;
    kmer_bar_ = rhs.kmer_bar_;
    kmer_hat_ptr = (rhs.kmer_hat_ptr == &rhs.kmer_ ? &kmer_ : &kmer_bar_);

    return *this;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::in_canonical_form() const
{
    return &kmer_ == kmer_hat_ptr;
}


template <uint16_t k>
inline void Directed_Vertex<k>::from_super_kmer(const uint64_t* const super_kmer, const std::size_t word_count)
{
    kmer_.from_super_kmer(super_kmer, word_count);
    init();
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::kmer() const
{
    return kmer_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::kmer_bar() const
{
    return kmer_bar_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::canonical() const
{
    return *kmer_hat_ptr;
}


template <uint16_t k>
inline void Directed_Vertex<k>::roll_forward(const DNA::Base b)
{
    kmer_.roll_to_next_kmer(b, kmer_bar_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
}


template <uint16_t k>
inline const Directed_Vertex<k> Directed_Vertex<k>::roll_forward(const DNA::Base b) const
{
    Directed_Vertex<k> temp(*this);

    temp.roll_forward(b);
    return temp;
}


template <uint16_t k>
inline void Directed_Vertex<k>::roll_backward(const DNA::Base b)
{
    kmer_.roll_to_prev_kmer(b, kmer_bar_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
}


template <uint16_t k>
inline const Directed_Vertex<k> Directed_Vertex<k>::roll_backward(const DNA::Base b) const
{
    Directed_Vertex temp(*this);

    temp.roll_backward(b);
    return temp;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::is_same_vertex(const Directed_Vertex<k>& v) const
{
    // TODO: revert back w/ correct design.
    // return hash() == v.hash();
    return canonical() == v.canonical();
}


}



#endif
