// Author: Ana Angelov
#include <cassert>
#include <cstdio>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "bamboo_filter.h"
#include "hash.h"

static std::string rand_dna(std::mt19937_64& rng, size_t len) {
    static const char* abc = "ACGT";
    std::string s; s.resize(len);
    for (size_t i = 0; i < len; ++i) s[i] = abc[rng() & 3];
    return s;
}

int main() {
    using namespace bamboo;

    {
        BambooConfig cfg;
        BambooFilter bf(cfg);
        std::mt19937_64 rng(13);
        std::vector<std::string> items;
        const int N = 200000;
        items.reserve(N);
        for (int i = 0; i < N; ++i) items.push_back(rand_dna(rng, 31));
        for (const auto& it : items) bf.insert(it);
        for (const auto& it : items) assert(bf.lookup(it));
        for (size_t i = 0; i < items.size() / 2; ++i) {
            assert(bf.remove(items[i]));
        }
        for (size_t i = items.size() / 2; i < items.size(); ++i) {
            assert(bf.lookup(items[i]));
        }
    }

    {
        BambooConfig cfg;
        cfg.bucket_per_segment = 64;
        cfg.initial_segments = 4;
        BambooFilter bf(cfg);
        for (int round = 0; round < 5; ++round) {
            bf.extend();
        }
        assert(bf.segment_count() > 4);
    }

    std::puts("test_bamboo OK");
    return 0;
}
