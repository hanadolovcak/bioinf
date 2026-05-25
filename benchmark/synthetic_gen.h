// Author: Hana Dolovčak
#ifndef BAMBOO_SYNTHETIC_GEN_H
#define BAMBOO_SYNTHETIC_GEN_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace bamboo {

// Returns a random DNA string of the given length using the given seed.
// The same seed always produces the same string.
std::string generate_random_dna(size_t length, uint32_t seed);

}

#endif
