// Author: Ana Angelov
#ifndef BAMBOO_HASH_H
#define BAMBOO_HASH_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace bamboo {

// Returns a hash value of a DNA k-mer.
// Uses ntHash. The k-mer must contain only DNA chars (A, C, G, T).
uint64_t hash_kmer(const std::string& kmer);

// Makes a fingerprint from the hash by extracting `fp_bits` bits
// starting at bit position `shift`.
uint32_t fingerprint_bits(uint64_t h, unsigned shift, unsigned fp_bits);

// Helper for finding another possible bucket from the given fingerprint.
uint64_t mini_hash(uint32_t fp);

class RollingKmerHasher {
public:
    // Sets up the hasher to walk over `sequence` with a window of length k.
    // The sequence string must live as long as the hasher.
    RollingKmerHasher(const std::string& sequence, size_t k);

    // Returns true if the window is currently on a valid k-mer.
    bool valid() const;

    // Moves the window one step to the right. Returns false if at the end.
    bool roll();

    // Returns the hash of the k-mer under the current window.
    uint64_t current_hash() const;

private:
    const char* seq_;
    size_t seq_len_;
    size_t k_;
    size_t pos_;
    uint64_t fh_val_;
    bool primed_;
};

}

#endif
