// Author: Hana Dolovčak
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

#include "kmer_extractor.h"

int main() {
    using namespace bamboo;

    {
        auto kmers = extract_kmers("ACGTACGT", 4);
        std::vector<std::string> expected = {"ACGT", "CGTA", "GTAC", "TACG", "ACGT"};
        assert(kmers == expected);
    }
    {
        auto kmers = extract_kmers("ACGT", 4);
        assert(kmers.size() == 1u);
        assert(kmers[0] == "ACGT");
    }
    {
        auto kmers = extract_kmers("ACG", 4);
        assert(kmers.empty());
    }
    {
        auto kmers = extract_kmers("ACGT", 0);
        assert(kmers.empty());
    }
    {
        size_t count = 0;
        for_each_kmer(std::string("ACGTACGT"), 4, [&](const char*, size_t) { ++count; });
        assert(count == 5u);
    }

    std::puts("test_kmer OK");
    return 0;
}
