// Author: Hana Dolovčak
#include "segment.h"

#include <cassert>

#include "hash.h"

namespace bamboo {

namespace {
constexpr int kMaxKicks = 500;
}

// Creates an empty segment with num_buckets x bucket_size slots.
// num_buckets must be a power of 2, fp_bits must be between 1 and 15.
// rng_seed seeds the RNG used for random choices during cuckoo eviction.
Segment::Segment(size_t num_buckets, size_t bucket_size, unsigned fp_bits,
                 uint64_t rng_seed)
    : num_buckets_(num_buckets),
      bucket_size_(bucket_size),
      fp_mask_(static_cast<uint16_t>((1u << fp_bits) - 1)),
      bucket_mask_(num_buckets - 1),
      data_(num_buckets * bucket_size, 0),
      occupancy_(0),
      rng_state_(rng_seed == 0 ? 0x9E3779B97F4A7C15ULL : rng_seed) {
    assert(num_buckets > 0 && (num_buckets & (num_buckets - 1)) == 0);
    assert(bucket_size > 0);
    assert(fp_bits >= 1 && fp_bits <= 15);
}

uint16_t Segment::read_slot(size_t bucket, size_t slot) const {
    return data_[bucket * bucket_size_ + slot];
}

void Segment::write_slot(size_t bucket, size_t slot, uint16_t value) {
    data_[bucket * bucket_size_ + slot] = value;
}

// Returns the alternative bucket index for fingerprint fp given the
// current bucket. Implements the partial-key cuckoo trick: h2 = h1 XOR
// mini_hash(fp), so either bucket can find the other from fp alone.
size_t Segment::alt_bucket(size_t bucket, uint32_t fp) const {
    return bucket ^ (mini_hash(fp) & bucket_mask_);
}

uint32_t Segment::next_random() const {
    uint64_t x = rng_state_;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng_state_ = x;
    return static_cast<uint32_t>((x * 0x2545F4914F6CDD1DULL) >> 32);
}

bool Segment::try_place(size_t bucket, uint32_t fp) {
    uint16_t encoded = static_cast<uint16_t>(kOccupiedBit | (fp & fp_mask_));
    for (size_t s = 0; s < bucket_size_; ++s) {
        if ((read_slot(bucket, s) & kOccupiedBit) == 0) {
            write_slot(bucket, s, encoded);
            ++occupancy_;
            return true;
        }
    }
    return false;
}

int Segment::find_in_bucket(size_t bucket, uint32_t fp) const {
    uint16_t target = static_cast<uint16_t>(kOccupiedBit | (fp & fp_mask_));
    for (size_t s = 0; s < bucket_size_; ++s) {
        if (read_slot(bucket, s) == target) {
            return static_cast<int>(s);
        }
    }
    return -1;
}

// Inserts fingerprint fp using h1 as the primary bucket index.
// Falls back to the alternative bucket, then to cuckoo eviction, then to
// the overflow list if eviction runs out of kicks.
// Returns true if placed in the main table, false if placed in overflow.
bool Segment::insert(uint32_t fp, size_t h1) {
    if (try_place(h1, fp)) return true;
    size_t h2 = alt_bucket(h1, fp);
    if (try_place(h2, fp)) return true;

    size_t cur_bucket = (next_random() & 1u) ? h2 : h1;
    uint32_t cur_fp = fp;
    for (int kicks = 0; kicks < kMaxKicks; ++kicks) {
        size_t slot = next_random() % bucket_size_;
        uint16_t evicted = read_slot(cur_bucket, slot);
        uint16_t encoded = static_cast<uint16_t>(kOccupiedBit | (cur_fp & fp_mask_));
        write_slot(cur_bucket, slot, encoded);
        cur_fp = static_cast<uint32_t>(evicted & fp_mask_);
        cur_bucket = alt_bucket(cur_bucket, cur_fp);
        if (try_place(cur_bucket, cur_fp)) {
            return true;
        }
    }
    // Eviction limit hit: park the orphan in the overflow list so it is not lost.
    overflow_.push_back(Entry{cur_bucket, cur_fp});
    ++occupancy_;
    return false;
}

// Returns true if fingerprint fp is present at bucket h1 or its alternative
// (or in the overflow list).
bool Segment::lookup(uint32_t fp, size_t h1) const {
    if (find_in_bucket(h1, fp) >= 0) return true;
    size_t h2 = alt_bucket(h1, fp);
    if (find_in_bucket(h2, fp) >= 0) return true;
    if (!overflow_.empty()) {
        uint32_t fp_value = fp & fp_mask_;
        for (const Entry& e : overflow_) {
            if (e.fp == fp_value && (e.bucket == h1 || e.bucket == h2)) {
                return true;
            }
        }
    }
    return false;
}

// Removes one occurrence of fingerprint fp at bucket h1 or its alternative
// (or from the overflow list). Returns true if something was removed.
bool Segment::remove(uint32_t fp, size_t h1) {
    int slot = find_in_bucket(h1, fp);
    if (slot >= 0) {
        write_slot(h1, static_cast<size_t>(slot), 0);
        --occupancy_;
        return true;
    }
    size_t h2 = alt_bucket(h1, fp);
    slot = find_in_bucket(h2, fp);
    if (slot >= 0) {
        write_slot(h2, static_cast<size_t>(slot), 0);
        --occupancy_;
        return true;
    }
    if (!overflow_.empty()) {
        uint32_t fp_value = fp & fp_mask_;
        for (size_t i = 0; i < overflow_.size(); ++i) {
            const Entry& e = overflow_[i];
            if (e.fp == fp_value && (e.bucket == h1 || e.bucket == h2)) {
                overflow_.erase(overflow_.begin() + static_cast<long>(i));
                --occupancy_;
                return true;
            }
        }
    }
    return false;
}

// occupancy() / capacity(), between 0 and 1.
double Segment::load_factor() const {
    return static_cast<double>(occupancy_) / static_cast<double>(capacity());
}

// Rough estimate of memory used by this segment in bytes.
size_t Segment::memory_bytes() const {
    return sizeof(Segment)
         + data_.capacity() * sizeof(uint16_t)
         + overflow_.capacity() * sizeof(Entry);
}

// Removes and returns all main-table entries whose fingerprint has
// bit `bit_index` equal to `bit_value`. Used by BambooFilter::extend
// to move half of a split segment's entries into the new segment.
std::vector<Segment::Entry> Segment::drain_entries_with_fp_bit(
    unsigned bit_index, bool bit_value) {
    std::vector<Entry> drained;
    uint16_t mask = static_cast<uint16_t>(1u << bit_index);
    for (size_t b = 0; b < num_buckets_; ++b) {
        for (size_t s = 0; s < bucket_size_; ++s) {
            uint16_t v = data_[b * bucket_size_ + s];
            if ((v & kOccupiedBit) == 0) continue;
            uint32_t fp = static_cast<uint32_t>(v & fp_mask_);
            bool bit_set = (v & mask) != 0;
            if (bit_set == bit_value) {
                drained.push_back({b, fp});
                data_[b * bucket_size_ + s] = 0;
                --occupancy_;
            }
        }
    }
    return drained;
}

// Same as drain_entries_with_fp_bit but for the overflow list.
std::vector<Segment::Entry> Segment::drain_overflow_with_fp_bit(
    unsigned bit_index, bool bit_value) {
    std::vector<Entry> drained;
    uint32_t mask = 1u << bit_index;
    auto it = overflow_.begin();
    while (it != overflow_.end()) {
        bool bit_set = (it->fp & mask) != 0;
        if (bit_set == bit_value) {
            drained.push_back(*it);
            it = overflow_.erase(it);
            --occupancy_;
        } else {
            ++it;
        }
    }
    return drained;
}

}
