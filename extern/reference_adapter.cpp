// Author: Ana Angelov
//
// Adapter that wraps the reference Bamboo Filter implementation
//     https://github.com/wanghanchengchn/bamboofilters
// inside the same CLI / CSV harness as our `benchmark` executable.
//
// Build prerequisite: the reference repository must be cloned next to this
// project at $REPO_ROOT/../bamboofilters (sibling of the bioinf project root).
// The reference uses x86 AVX2 intrinsics, so this adapter only builds on
// x86_64 with `-mavx2`. CMake gates the target accordingly.
//
// Each row this binary appends matches the schema produced by
// `benchmark_main.cpp` so that scripts/compare.py can diff them directly.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "bamboofilter/bamboofilter.hpp"
#include "bamboofilter/bitsutil.h"

#include "fasta_parser.h"
#include "memory_util.h"
#include "synthetic_gen.h"
#include "timer.h"

namespace {

struct Args {
    std::string mode = "synthetic";
    size_t k = 20;
    size_t data_size = 100000;
    std::string fasta_path;
    std::string output_path;
    uint32_t seed = 42;
    uint32_t initial_capacity = 0;  // 0 → upperpower2(num_unique)
    uint32_t split_bits = 2;
    size_t fpr_queries = 0;
};

void print_usage() {
    std::cerr <<
        "Usage:\n"
        "  reference_adapter --mode synthetic|ecoli [options] --output FILE.csv\n";
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << name << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (s == "--mode") a.mode = need("--mode");
        else if (s == "--k") a.k = std::strtoull(need("--k"), nullptr, 10);
        else if (s == "--data-size") a.data_size = std::strtoull(need("--data-size"), nullptr, 10);
        else if (s == "--fasta") a.fasta_path = need("--fasta");
        else if (s == "--seed") a.seed = static_cast<uint32_t>(std::strtoul(need("--seed"), nullptr, 10));
        else if (s == "--initial-capacity") a.initial_capacity = static_cast<uint32_t>(std::strtoul(need("--initial-capacity"), nullptr, 10));
        else if (s == "--split-bits") a.split_bits = static_cast<uint32_t>(std::strtoul(need("--split-bits"), nullptr, 10));
        else if (s == "--fpr-queries") a.fpr_queries = std::strtoull(need("--fpr-queries"), nullptr, 10);
        else if (s == "--output") a.output_path = need("--output");
        else if (s == "--help" || s == "-h") { print_usage(); std::exit(0); }
        else if (s == "--filter") { (void)need("--filter"); }  // accepted+ignored for symmetry
        else if (s == "--fp-bits") { (void)need("--fp-bits"); }  // hardcoded to 12 in their impl
        else if (s == "--initial-segments") { (void)need("--initial-segments"); }
        else if (s == "--bucket-per-segment") { (void)need("--bucket-per-segment"); }
        else { std::cerr << "unknown arg: " << s << "\n"; return false; }
    }
    return true;
}

std::string load_sequence(const Args& a) {
    if (a.mode == "synthetic") {
        return bamboo::generate_random_dna(a.data_size, a.seed);
    }
    if (a.mode == "ecoli") {
        if (a.fasta_path.empty()) {
            std::cerr << "--fasta required for ecoli mode\n";
            std::exit(2);
        }
        return bamboo::load_fasta(a.fasta_path);
    }
    std::cerr << "unknown mode: " << a.mode << "\n";
    std::exit(2);
}

// Collects all length-k substrings of `seq` into `out`, one std::string per
// k-mer. The reference Insert/Lookup take `const char*` so we need owning
// storage anyway.
void collect_kmers(const std::string& seq, size_t k,
                   std::vector<std::string>& out) {
    if (k == 0 || seq.size() < k) return;
    out.reserve(seq.size() - k + 1);
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        out.emplace_back(seq, i, k);
    }
}

}

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) return 2;

    std::string seq = load_sequence(args);
    if (seq.size() < args.k) {
        std::cerr << "sequence shorter than k\n";
        return 1;
    }

    std::vector<std::string> kmers;
    collect_kmers(seq, args.k, kmers);
    size_t num_kmers = kmers.size();

    std::unordered_set<std::string> seen;
    seen.reserve(num_kmers);
    std::vector<const std::string*> unique;
    unique.reserve(num_kmers);
    for (const auto& s : kmers) {
        if (seen.insert(s).second) unique.push_back(&s);
    }
    size_t num_unique = unique.size();

    uint32_t cap = args.initial_capacity;
    if (cap == 0) {
        cap = static_cast<uint32_t>(num_unique == 0 ? 1 : num_unique);
        cap = static_cast<uint32_t>(upperpower2(cap));
    }
    // The reference impl computes INIT_TABLE_BITS = ceil(log2(capacity/4)) and
    // then derives NUM_SEG_BITS = INIT_TABLE_BITS - BUCKETS_PER_SEG. With
    // BUCKETS_PER_SEG = 10 and capacity < 4096 the subtraction underflows
    // (uint32_t) and the constructor tries to allocate ~4 billion segments.
    // Clamp upward to avoid that.
    const uint32_t kMinCap = 4u << BUCKETS_PER_SEG;
    if (cap < kMinCap) cap = kMinCap;

    BambooFilter* bbf = new BambooFilter(cap, args.split_bits);

    bamboo::Timer t_insert;
    for (const std::string* sp : unique) {
        bbf->Insert(sp->c_str());
    }
    double insert_s = t_insert.seconds();

    size_t peak_kb = bamboo::peak_rss_kb();

    bamboo::Timer t_lookup;
    size_t pos_hits = 0;
    for (const std::string* sp : unique) {
        if (bbf->Lookup(sp->c_str())) ++pos_hits;
    }
    double lookup_s = t_lookup.seconds();

    double fpr = -1.0;
    if (args.fpr_queries > 0) {
        std::mt19937_64 rng(static_cast<uint64_t>(args.seed) ^ 0xA5A5A5A5ULL);
        size_t false_positives = 0;
        size_t trials = 0;
        std::string buf(args.k, 'A');
        const char abc[4] = {'A', 'C', 'G', 'T'};
        while (trials < args.fpr_queries) {
            for (size_t i = 0; i < args.k; ++i) buf[i] = abc[rng() & 3u];
            if (seen.count(buf)) continue;
            ++trials;
            if (bbf->Lookup(buf.c_str())) ++false_positives;
        }
        fpr = static_cast<double>(false_positives) / static_cast<double>(trials);
    }

    // The reference implementation does not expose its internal byte count
    // directly. We approximate using peak_rss minus an empirically-small
    // baseline (~5 MB on Linux). When in doubt, compare on peak_rss only.
    size_t internal_bytes = 0;

    std::ostringstream row;
    row << args.mode << ','
        << (args.mode == "ecoli" ? args.fasta_path : "synthetic") << ','
        << seq.size() << ','
        << args.k << ','
        << num_kmers << ','
        << num_unique << ','
        << "reference" << ','
        << insert_s << ','
        << lookup_s << ','
        << fpr << ','
        << peak_kb << ','
        << internal_bytes << ','
        << bbf->hash_table_.size() << ','
        << 0.0 << ','  // load factor — reference does not expose
        << pos_hits;

    if (!args.output_path.empty()) {
        bool need_header = false;
        {
            std::ifstream check(args.output_path);
            need_header = !check.good();
        }
        std::ofstream out(args.output_path, std::ios::app);
        if (need_header) {
            out << "mode,data_source,data_size,k,num_kmers,num_unique,filter,"
                   "insert_s,lookup_s,fpr,peak_rss_kb,internal_bytes,"
                   "segments,load_factor,positive_hits\n";
        }
        out << row.str() << "\n";
    }

    std::cout << "REFERENCE mode=" << args.mode
              << " k=" << args.k
              << " kmers=" << num_kmers
              << " unique=" << num_unique
              << " insert_s=" << insert_s
              << " lookup_s=" << lookup_s
              << " fpr=" << fpr
              << " peak_kb=" << peak_kb
              << " segs=" << bbf->hash_table_.size()
              << " pos_hits=" << pos_hits
              << "\n";

    delete bbf;
    return 0;
}
