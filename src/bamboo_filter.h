// Author: Ana Angelov
#ifndef BAMBOO_BAMBOO_FILTER_H
#define BAMBOO_BAMBOO_FILTER_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "segment.h"

namespace bamboo {

// Configuration for BambooFilter. All thresholds are between 0 and 1.
struct BambooConfig {
    // Initial number of segments. Must be a power of 2.
    size_t initial_segments    = 16;
    // Buckets per segment. Must be a power of 2.
    size_t bucket_per_segment  = 256;
    // Slots per bucket.
    size_t bucket_size         = 4;
    // Fingerprint bit width (between 1 and 15).
    unsigned fp_bits           = 8;
    // Average load that triggers automatic extend().
    double extend_threshold    = 0.95;
    // Average load below which automatic compress() runs after a remove.
    double compress_threshold  = 0.40;
    // Disables automatic extend (useful for tests).
    bool enable_extend         = true;
    // Disables automatic compress (useful for tests).
    bool enable_compress       = true;
};

class BambooFilter {
public:
    // Builds a Bamboo filter using config. Default config gives 16 segments
    // of 256 x 4 slots with 8-bit fingerprints.
    explicit BambooFilter(const BambooConfig& config = {});

    // Inserts a DNA k-mer into the filter.
    void insert(const std::string& item);

    // Returns true if item is (probably) in the filter. May give false
    // positives at a rate controlled by fp_bits.
    bool lookup(const std::string& item) const;

    // Removes one occurrence of item. Returns true if something was removed.
    // May trigger a compress() if the load drops below the compress threshold.
    bool remove(const std::string& item);

    // Same as insert but takes a pre-computed hash. Useful when the caller
    // already has the hash (e.g. from RollingKmerHasher or in benchmarks).
    void insert_hash(uint64_t h);

    // Same as lookup but takes a pre-computed hash.
    bool lookup_hash(uint64_t h) const;

    // Number of items inserted (counts duplicates).
    size_t size() const { return total_inserted_; }

    // Number of segments currently allocated.
    size_t segment_count() const { return segments_.size(); }

    // Overall load factor across all segments.
    double load_factor() const;

    // Rough estimate of memory used by the filter in bytes.
    size_t memory_bytes() const;

    // Splits the next segment in two and moves half of its entries to a
    // newly allocated segment. Normally called automatically by insert_hash;
    // exposed publicly so tests can force a split.
    void extend();

private:
    // Computes (segment index, bucket within segment, fingerprint) from a hash.
    void decompose(uint64_t h, size_t& seg, size_t& bucket, uint32_t& fp) const;

    // Inverse of extend(): merges the last segment back into the one that
    // was split last. Auto-triggered by remove() when load drops below
    // config_.compress_threshold.
    void compress();

    BambooConfig config_;
    unsigned l_b_;
    unsigned l_s0_;
    unsigned expand_round_;
    size_t split_pointer_;
    size_t n_initial_segments_;
    size_t total_inserted_;
    size_t expand_trigger_;
    uint64_t rng_seed_counter_;

    std::vector<std::unique_ptr<Segment>> segments_;
};

}

#endif
