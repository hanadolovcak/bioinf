// Author: Hana Dolovčak
//
// CLI entry point for the Bamboo filter benchmark.
// Reads a sequence (random DNA or a FASTA file), extracts k-mers, deduplicates
// them, inserts them into the filter and times the insert and lookup phases.
// Optionally estimates the false positive rate. Writes a CSV row with the
// measured numbers.
//
// Run with --help for the list of options.

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>

#include "bamboo_filter.h"
#include "fasta_parser.h"
#include "hash.h"
#include "memory_util.h"
#include "synthetic_gen.h"
#include "timer.h"

namespace {

// All CLI options parsed from argv.
struct Args {
    std::string mode = "synthetic";
    size_t k = 20;
    size_t data_size = 100000;
    std::string filter = "bamboo";
    std::string fasta_path;
    std::string output_path;
    uint32_t seed = 42;
    size_t fp_bits = 8;
    size_t initial_segments = 16;
    size_t bucket_per_segment = 256;
    size_t fpr_queries = 0;
};

// Prints the CLI usage text to stderr.
void print_usage() {
    std::cerr <<
        "Usage:\n"
        "  benchmark --mode synthetic|ecoli [options] --output FILE.csv\n"
        "\n"
        "Options:\n"
        "  --mode {synthetic|ecoli}            workload type (default synthetic)\n"
        "  --k N                                k-mer length (default 20)\n"
        "  --data-size N                        synthetic sequence length (default 100000)\n"
        "  --fasta PATH                         FASTA path (required for ecoli mode)\n"
        "  --filter {bamboo}                    filter under test (default bamboo)\n"
        "  --seed S                             RNG seed (default 42)\n"
        "  --fp-bits N                          fingerprint bits 1..15 (default 8)\n"
        "  --initial-segments N                 initial bamboo segment count, power of two\n"
        "  --bucket-per-segment N               buckets per segment, power of two\n"
        "  --fpr-queries N                      negative queries for FPR estimation\n"
        "  --output PATH                        CSV file to append results to\n";
}

// Fills `a` from argv. Returns false on an unknown flag.
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
        else if (s == "--filter") a.filter = need("--filter");
        else if (s == "--seed") a.seed = static_cast<uint32_t>(std::strtoul(need("--seed"), nullptr, 10));
        else if (s == "--fp-bits") a.fp_bits = std::strtoull(need("--fp-bits"), nullptr, 10);
        else if (s == "--initial-segments") a.initial_segments = std::strtoull(need("--initial-segments"), nullptr, 10);
        else if (s == "--bucket-per-segment") a.bucket_per_segment = std::strtoull(need("--bucket-per-segment"), nullptr, 10);
        else if (s == "--fpr-queries") a.fpr_queries = std::strtoull(need("--fpr-queries"), nullptr, 10);
        else if (s == "--output") a.output_path = need("--output");
        else if (s == "--help" || s == "-h") { print_usage(); std::exit(0); }
        else { std::cerr << "unknown arg: " << s << "\n"; return false; }
    }
    return true;
}

// Returns the DNA sequence to process based on the chosen mode.
// For synthetic mode it generates a random string; for ecoli mode it loads
// the FASTA at a.fasta_path.
std::string load_sequence(const Args& a) {
    if (a.mode == "synthetic") {
        return bamboo::generate_random_dna(a.data_size, a.seed);
    } else if (a.mode == "ecoli") {
        if (a.fasta_path.empty()) {
            std::cerr << "--fasta required for ecoli mode\n";
            std::exit(2);
        }
        return bamboo::load_fasta(a.fasta_path);
    } else {
        std::cerr << "unknown mode: " << a.mode << "\n";
        std::exit(2);
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

    size_t num_kmers = seq.size() - args.k + 1;

    bamboo::BambooConfig config;
    config.bucket_per_segment = args.bucket_per_segment;
    config.initial_segments = args.initial_segments;
    config.fp_bits = static_cast<unsigned>(args.fp_bits);
    bamboo::BambooFilter bf(config);

    std::unordered_set<uint64_t> unique_hashes;
    unique_hashes.reserve(num_kmers);
    {
        bamboo::RollingKmerHasher hasher(seq, args.k);
        if (hasher.valid()) {
            do {
                unique_hashes.insert(hasher.current_hash());
            } while (hasher.roll());
        }
    }
    size_t num_unique = unique_hashes.size();

    bamboo::Timer t_insert;
    for (uint64_t h : unique_hashes) {
        bf.insert_hash(h);
    }
    double insert_s = t_insert.seconds();

    size_t peak_kb = bamboo::peak_rss_kb();

    bamboo::Timer t_lookup;
    size_t pos_hits = 0;
    for (uint64_t h : unique_hashes) {
        if (bf.lookup_hash(h)) ++pos_hits;
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
            uint64_t h = bamboo::hash_kmer(buf);
            if (unique_hashes.count(h)) continue;
            ++trials;
            if (bf.lookup_hash(h)) ++false_positives;
        }
        fpr = static_cast<double>(false_positives) / static_cast<double>(trials);
    }

    size_t internal_bytes = bf.memory_bytes();

    std::ostringstream row;
    row << args.mode << ','
        << (args.mode == "ecoli" ? args.fasta_path : "synthetic") << ','
        << seq.size() << ','
        << args.k << ','
        << num_kmers << ','
        << num_unique << ','
        << args.filter << ','
        << insert_s << ','
        << lookup_s << ','
        << fpr << ','
        << peak_kb << ','
        << internal_bytes << ','
        << bf.segment_count() << ','
        << bf.load_factor() << ','
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

    std::cout << "mode=" << args.mode
              << " k=" << args.k
              << " kmers=" << num_kmers
              << " unique=" << num_unique
              << " insert_s=" << insert_s
              << " lookup_s=" << lookup_s
              << " fpr=" << fpr
              << " peak_kb=" << peak_kb
              << " internal_bytes=" << internal_bytes
              << " segs=" << bf.segment_count()
              << " load=" << bf.load_factor()
              << " pos_hits=" << pos_hits
              << "\n";

    return 0;
}
