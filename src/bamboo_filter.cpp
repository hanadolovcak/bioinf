// Author: Ana Angelov
#include "bamboo_filter.h"

#include <cassert>

#include "hash.h"

namespace bamboo {

namespace {

bool is_power_of_two(size_t x) {
    return x > 0 && (x & (x - 1)) == 0;
}

unsigned log2_exact(size_t x) {
    unsigned r = 0;
    while (x > 1) { x >>= 1; ++r; }
    return r;
}

}

// Builds a Bamboo filter using config. Default config gives 16 segments of
// 256 x 4 slots with 8-bit fingerprints.
BambooFilter::BambooFilter(const BambooConfig& config)
    : config_(config),
      expand_round_(0),
      split_pointer_(0),
      n_initial_segments_(config.initial_segments),
      total_inserted_(0),
      rng_seed_counter_(0xDEADBEEFULL) {
    assert(is_power_of_two(config.bucket_per_segment));
    assert(is_power_of_two(config.initial_segments));
    assert(config.bucket_size > 0);
    assert(config.fp_bits >= 1 && config.fp_bits <= 15);

    l_b_ = log2_exact(config.bucket_per_segment);
    l_s0_ = log2_exact(config.initial_segments);

    double k = static_cast<double>(config.bucket_per_segment)
             * static_cast<double>(config.bucket_size)
             * config.extend_threshold;
    expand_trigger_ = (k < 1.0) ? 1 : static_cast<size_t>(k);

    segments_.reserve(config.initial_segments * 2);
    for (size_t i = 0; i < config.initial_segments; ++i) {
        segments_.push_back(std::make_unique<Segment>(
            config.bucket_per_segment, config.bucket_size, config.fp_bits,
            rng_seed_counter_++));
    }
}

// Computes (segment index, bucket within segment, fingerprint) from a hash.
void BambooFilter::decompose(uint64_t h, size_t& seg, size_t& bucket,
                             uint32_t& fp) const {
    bucket = static_cast<size_t>(h & ((1ULL << l_b_) - 1));
    fp = fingerprint_bits(h, l_b_ + l_s0_, config_.fp_bits);

    size_t small_count = (static_cast<size_t>(1) << expand_round_)
                       * n_initial_segments_;
    size_t large_count = small_count << 1;
    size_t cand = static_cast<size_t>((h >> l_b_) & (large_count - 1));

    if (cand < small_count) {
        seg = cand;
    } else {
        size_t offset = cand - small_count;
        seg = (offset < split_pointer_) ? cand : offset;
    }
}

// Same as insert but takes a pre-computed hash. Useful when the caller
// already has the hash (e.g. from RollingKmerHasher or in benchmarks).
void BambooFilter::insert_hash(uint64_t h) {
    if (config_.enable_extend && expand_trigger_ > 0
        && total_inserted_ > 0
        && (total_inserted_ % expand_trigger_) == 0) {
        extend();
    }

    size_t seg, bucket;
    uint32_t fp;
    decompose(h, seg, bucket, fp);

    if (!segments_[seg]->insert(fp, bucket)) {
        // Insert overflowed — trigger an extend to spread the load.
        if (config_.enable_extend) {
            extend();
        }
    }
    ++total_inserted_;
}

// Same as lookup but takes a pre-computed hash.
bool BambooFilter::lookup_hash(uint64_t h) const {
    size_t seg, bucket;
    uint32_t fp;
    decompose(h, seg, bucket, fp);
    return segments_[seg]->lookup(fp, bucket);
}

// Inserts a DNA k-mer into the filter.
void BambooFilter::insert(const std::string& item) {
    insert_hash(hash_kmer(item));
}

// Returns true if item is (probably) in the filter. May give false
// positives at a rate controlled by fp_bits.
bool BambooFilter::lookup(const std::string& item) const {
    return lookup_hash(hash_kmer(item));
}

// Removes one occurrence of item. Returns true if something was removed.
// May trigger a compress() if the load drops below the compress threshold.
bool BambooFilter::remove(const std::string& item) {
    uint64_t h = hash_kmer(item);
    size_t seg, bucket;
    uint32_t fp;
    decompose(h, seg, bucket, fp);
    if (!segments_[seg]->remove(fp, bucket)) {
        return false;
    }
    if (total_inserted_ > 0) --total_inserted_;

    if (config_.enable_compress
        && segments_.size() > n_initial_segments_
        && load_factor() < config_.compress_threshold) {
        compress();
    }
    return true;
}

// Splits the next segment in two and moves half of its entries to a
// newly allocated segment. Normally called automatically by insert_hash;
// exposed publicly so tests can force a split.
void BambooFilter::extend() {
    if (expand_round_ >= config_.fp_bits) {
        return;  // fingerprint bits exhausted; cannot split deterministically further
    }

    size_t p = split_pointer_;
    size_t new_seg_idx = segments_.size();

    auto new_seg = std::make_unique<Segment>(
        config_.bucket_per_segment, config_.bucket_size, config_.fp_bits,
        rng_seed_counter_++);

    auto drained = segments_[p]->drain_entries_with_fp_bit(
        expand_round_, /*bit_value=*/true);
    auto drained_overflow = segments_[p]->drain_overflow_with_fp_bit(
        expand_round_, /*bit_value=*/true);

    segments_.push_back(std::move(new_seg));

    for (const auto& e : drained) {
        segments_[new_seg_idx]->insert(e.fp, e.bucket);
    }
    for (const auto& e : drained_overflow) {
        segments_[new_seg_idx]->insert(e.fp, e.bucket);
    }

    size_t small_count = (static_cast<size_t>(1) << expand_round_)
                       * n_initial_segments_;
    ++split_pointer_;
    if (split_pointer_ >= small_count) {
        split_pointer_ = 0;
        ++expand_round_;
    }
}

// Inverse of extend(): merges the last segment back into the one that
// was split last. Auto-triggered by remove() when load drops below
// config_.compress_threshold.
void BambooFilter::compress() {
    if (segments_.size() <= n_initial_segments_) return;

    if (split_pointer_ == 0) {
        if (expand_round_ == 0) return;  // already at minimum size
        --expand_round_;
        size_t small_count = (static_cast<size_t>(1) << expand_round_)
                           * n_initial_segments_;
        split_pointer_ = small_count;
    }
    --split_pointer_;

    size_t p = split_pointer_;
    size_t last_idx = segments_.size() - 1;

    segments_[last_idx]->for_each_entry([&](Segment::Entry e) {
        segments_[p]->insert(e.fp, e.bucket);
    });

    segments_.pop_back();
}

// Overall load factor across all segments.
double BambooFilter::load_factor() const {
    if (segments_.empty()) return 0.0;
    size_t cap = 0, occ = 0;
    for (const auto& s : segments_) {
        cap += s->capacity();
        occ += s->occupancy();
    }
    return static_cast<double>(occ) / static_cast<double>(cap);
}

// Rough estimate of memory used by the filter in bytes.
size_t BambooFilter::memory_bytes() const {
    size_t total = sizeof(*this);
    for (const auto& s : segments_) {
        total += s->memory_bytes();
    }
    total += segments_.capacity() * sizeof(std::unique_ptr<Segment>);
    return total;
}

}
