// Author: Hana Dolovčak
#ifndef BAMBOO_KMER_EXTRACTOR_H
#define BAMBOO_KMER_EXTRACTOR_H

#include <cstddef>
#include <string>
#include <vector>

namespace bamboo {

// Returns all length-k contiguous substrings (k-mers) of seq.
// Returns an empty vector if k == 0 or k > seq.size().
std::vector<std::string> extract_kmers(const std::string& seq, size_t k);

// Calls callback(data, k) for each k-mer in seq, without allocating strings.
// Useful for very long sequences (e.g. whole genomes).
template <typename F>
void for_each_kmer(const std::string& seq, size_t k, F&& callback) {
    if (k == 0 || k > seq.size()) return;
    const char* data = seq.data();
    size_t last_start = seq.size() - k;
    for (size_t i = 0; i <= last_start; ++i) {
        callback(data + i, k);
    }
}

}

#endif
