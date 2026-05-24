// Author: Ana Angelov
#include "hash.h"

#include "nthash.hpp"

namespace bamboo {

// Returns a hash value of a DNA k-mer.
// Uses ntHash. The k-mer must contain only DNA chars (A, C, G, T).
uint64_t hash_kmer(const std::string& kmer) {
    return NTF64(kmer.data(), static_cast<unsigned>(kmer.size()));
}

// Makes a fingerprint from the hash by extracting `fp_bits` bits
// starting at bit position `shift`.
uint32_t fingerprint_bits(uint64_t h, unsigned shift, unsigned fp_bits) {
    return static_cast<uint32_t>((h >> shift) & ((1ULL << fp_bits) - 1));
}

// Helper for finding another possible bucket from the given fingerprint.
uint64_t mini_hash(uint32_t fp) {
    uint64_t z = static_cast<uint64_t>(fp) + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return z;
}

// Sets up the hasher to walk over `sequence` with a window of length k.
// The sequence string must live as long as the hasher.
RollingKmerHasher::RollingKmerHasher(const std::string& sequence, size_t k)
    : seq_(sequence.data()),
      seq_len_(sequence.size()),
      k_(k),
      pos_(0),
      fh_val_(0),
      primed_(false) {
    if (k_ > 0 && seq_len_ >= k_) {
        fh_val_ = NTF64(seq_, static_cast<unsigned>(k_));
        primed_ = true;
    }
}

// Returns true if the window is currently on a valid k-mer.
bool RollingKmerHasher::valid() const {
    return primed_ && pos_ + k_ <= seq_len_;
}

// Moves the window one step to the right. Returns false if at the end.
bool RollingKmerHasher::roll() {
    if (!primed_ || pos_ + k_ >= seq_len_) {
        return false;
    }
    unsigned char char_out = static_cast<unsigned char>(seq_[pos_]);
    unsigned char char_in = static_cast<unsigned char>(seq_[pos_ + k_]);
    ++pos_;
    fh_val_ = NTF64(fh_val_, static_cast<unsigned>(k_), char_out, char_in);
    return true;
}

// Returns the hash of the k-mer under the current window.
uint64_t RollingKmerHasher::current_hash() const {
    return fh_val_;
}

}
