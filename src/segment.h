// Author: Hana Dolovčak
#ifndef BAMBOO_SEGMENT_H
#define BAMBOO_SEGMENT_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace bamboo {

class Segment {
public:
    // A single stored entry: which bucket it lives in and its fingerprint.
    struct Entry {
        size_t bucket;
        uint32_t fp;
    };

    // Creates an empty segment with num_buckets x bucket_size slots.
    // num_buckets must be a power of 2, fp_bits must be between 1 and 15.
    // rng_seed seeds the RNG used for random choices during cuckoo eviction.
    Segment(size_t num_buckets, size_t bucket_size, unsigned fp_bits,
            uint64_t rng_seed);

    // Inserts fingerprint fp using h1 as the primary bucket index.
    // Falls back to the alternative bucket, then to cuckoo eviction, then to
    // the overflow list if eviction runs out of kicks.
    // Returns true if placed in the main table, false if placed in overflow.
    bool insert(uint32_t fp, size_t h1);

    // Returns true if fingerprint fp is present at bucket h1 or its alternative
    // (or in the overflow list).
    bool lookup(uint32_t fp, size_t h1) const;

    // Removes one occurrence of fingerprint fp at bucket h1 or its alternative
    // (or from the overflow list). Returns true if something was removed.
    bool remove(uint32_t fp, size_t h1);

    // Returns the alternative bucket index for fingerprint fp given the
    // current bucket. Implements the partial-key cuckoo trick: h2 = h1 XOR
    // mini_hash(fp), so either bucket can find the other from fp alone.
    size_t alt_bucket(size_t bucket, uint32_t fp) const;

    // Number of occupied slots (main table + overflow).
    size_t occupancy() const { return occupancy_; }

    // Total slots in the main table (excludes overflow).
    size_t capacity() const { return num_buckets_ * bucket_size_; }

    // occupancy() / capacity(), between 0 and 1.
    double load_factor() const;

    // Rough estimate of memory used by this segment in bytes.
    size_t memory_bytes() const;

    // Removes and returns all main-table entries whose fingerprint has
    // bit `bit_index` equal to `bit_value`. Used by BambooFilter::extend
    // to move half of a split segment's entries into the new segment.
    std::vector<Entry> drain_entries_with_fp_bit(unsigned bit_index, bool bit_value);

    // Same as drain_entries_with_fp_bit but for the overflow list.
    std::vector<Entry> drain_overflow_with_fp_bit(unsigned bit_index, bool bit_value);

    // Calls cb(entry) for every stored entry (main table and overflow).
    // Read-only; does not modify the segment.
    template <typename Callback>
    void for_each_entry(Callback&& cb) const;

private:
    static constexpr uint16_t kOccupiedBit = 0x8000;

    size_t num_buckets_;
    size_t bucket_size_;
    uint16_t fp_mask_;
    size_t bucket_mask_;
    std::vector<uint16_t> data_;
    std::vector<Entry> overflow_;
    size_t occupancy_;
    mutable uint64_t rng_state_;

    uint16_t read_slot(size_t bucket, size_t slot) const;
    void write_slot(size_t bucket, size_t slot, uint16_t value);
    bool try_place(size_t bucket, uint32_t fp);
    int find_in_bucket(size_t bucket, uint32_t fp) const;
    uint32_t next_random() const;
};

template <typename Callback>
void Segment::for_each_entry(Callback&& cb) const {
    for (size_t b = 0; b < num_buckets_; ++b) {
        for (size_t s = 0; s < bucket_size_; ++s) {
            uint16_t v = data_[b * bucket_size_ + s];
            if (v & kOccupiedBit) {
                cb(Entry{b, static_cast<uint32_t>(v & fp_mask_)});
            }
        }
    }
    for (const Entry& e : overflow_) {
        cb(e);
    }
}

}

#endif
