// Author: Hana Dolovčak
#include "synthetic_gen.h"

#include <random>

namespace bamboo {

// Returns a random DNA string of the given length using the given seed.
// The same seed always produces the same string.
std::string generate_random_dna(size_t length, uint32_t seed) {
    static const char alphabet[4] = {'A', 'C', 'G', 'T'};
    std::mt19937_64 rng(static_cast<uint64_t>(seed));
    std::string out;
    out.resize(length);
    for (size_t i = 0; i < length; ++i) {
        out[i] = alphabet[rng() & 3u];
    }
    return out;
}

}
