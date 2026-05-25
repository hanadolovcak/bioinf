// Author: Hana Dolovčak
#include <cassert>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "hash.h"
#include "segment.h"

static std::string rand_dna(std::mt19937_64& rng, size_t len) {
    static const char* abc = "ACGT";
    std::string s; s.resize(len);
    for (size_t i = 0; i < len; ++i) s[i] = abc[rng() & 3];
    return s;
}

int main() {
    using namespace bamboo;

    {
        Segment seg(256, 4, 8, 1);
        std::mt19937_64 rng(7);
        std::vector<std::pair<size_t, uint32_t>> inserted;
        const int N = 800;
        for (int i = 0; i < N; ++i) {
            uint64_t h = hash_kmer(rand_dna(rng, 24));
            size_t b = h & 0xFF;
            uint32_t fp = static_cast<uint32_t>((h >> 12) & 0xFF);
            assert(seg.insert(fp, b));
            inserted.push_back({b, fp});
        }
        for (auto p : inserted) {
            assert(seg.lookup(p.second, p.first));
        }
        for (size_t i = 0; i < inserted.size() / 2; ++i) {
            assert(seg.remove(inserted[i].second, inserted[i].first));
        }
        for (size_t i = inserted.size() / 2; i < inserted.size(); ++i) {
            assert(seg.lookup(inserted[i].second, inserted[i].first));
        }
    }

    {
        // fp=0 lifecycle
        Segment seg(64, 4, 8, 99);
        assert(seg.insert(0, 10));
        assert(seg.lookup(0, 10));
        assert(seg.remove(0, 10));
        assert(!seg.lookup(0, 10));
    }

    {
        // alt-bucket symmetry: inserting at h1, lookup via h2
        Segment seg(64, 4, 8, 99);
        uint32_t fp = 42;
        size_t b1 = 5;
        assert(seg.insert(fp, b1));
        size_t b2 = seg.alt_bucket(b1, fp);
        assert(b1 != b2);
        assert(seg.lookup(fp, b2));
    }

    std::puts("test_segment OK");
    return 0;
}
