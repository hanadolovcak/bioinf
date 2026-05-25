// Author: Hana Dolovčak
#include "kmer_extractor.h"

namespace bamboo {

// Returns all length-k contiguous substrings (k-mers) of seq.
// Returns an empty vector if k == 0 or k > seq.size().
std::vector<std::string> extract_kmers(const std::string& seq, size_t k) {
    std::vector<std::string> out;
    if (k == 0 || k > seq.size()) return out;
    out.reserve(seq.size() - k + 1);
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        out.emplace_back(seq, i, k);
    }
    return out;
}

}
