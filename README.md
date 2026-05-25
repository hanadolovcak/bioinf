# Bamboo Filter (C++17)

A from-scratch C++17 implementation of the **Bamboo Filter** approximate-membership
data structure (Wang et al., IEEE/ACM Transactions on Networking, 2024), tested on
synthetic DNA and on the *Escherichia coli* K-12 genome with k-mers of varying
length.

Project for the **Bioinformatika 1** course (FER, 2025/2026), task (2) — MDL
track. Course page: <https://www.fer.unizg.hr/predmet/bio1>.

## Authors

- Ana Angelov
- Hana Dolovčak

## Requirements

- Linux (tested on Ubuntu 22.04; should also work on Debian, Fedora, CentOS)
- CMake ≥ 3.16
- A C++17 compiler — g++ ≥ 9 or clang++ ≥ 10
- Python 3 with `matplotlib` (only needed for generating plots)


## Quickstart (Ubuntu / Debian)

```sh
# 1. Install build dependencies
sudo apt update
sudo apt install -y build-essential cmake git python3 python3-matplotlib

# 2. Clone and build
git clone https://github.com/hanadolovcak/bioinf.git bioinf
cd bioinf
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4

# 3. Run the test suite
ctest --test-dir build --output-on-failure

# 4. Run the full benchmark sweeps
chmod +x scripts/*.sh
scripts/run_synthetic.sh        # synthetic DNA, 10^3 to 10^7 bases, 5 k values
scripts/run_ecoli.sh            # E. coli K-12 genome, 5 k values

# 5. Generate figures from the results
python3 scripts/plot_results.py --csv results/synthetic.csv --prefix synthetic
python3 scripts/plot_results.py --csv results/ecoli.csv     --prefix ecoli
```

After this, `results/` contains both CSVs and `results/figures/` contains
the throughput, memory and FPR plots.

For a single benchmark run with custom parameters:

```sh
./build/benchmark --mode synthetic \
    --data-size 1000000 --k 20 \
    --fp-bits 12 --initial-segments 2048 \
    --fpr-queries 50000 \
    --output results/custom.csv
```

The header row of the CSV is written automatically the first time the file is
created; subsequent runs append one row each.

## CLI reference

```
--mode {synthetic|ecoli}      workload type (default synthetic)
--k N                         k-mer length (default 20)
--data-size N                 synthetic sequence length (default 100000)
--fasta PATH                  FASTA path (required for ecoli mode)
--seed S                      RNG seed (default 42)
--fp-bits N                   fingerprint bits 1..15 (default 8)
--initial-segments N          initial Bamboo segment count, must be a power of two
--bucket-per-segment N        buckets per segment, must be a power of two
--fpr-queries N               number of negative queries used to estimate FPR
--output PATH                 CSV file to append the row to
```

## Comparing with the reference implementation

The reference Bamboo Filter from
[wanghanchengchn/bamboofilters](https://github.com/wanghanchengchn/bamboofilters)
is a C++11 library without a matching CLI. We bundle a small adapter
(`extern/reference_adapter.cpp`) that wraps *their* `BambooFilter` class behind
the same CLI flags and CSV schema as our `benchmark`.

**Platform constraint.** The reference implementation uses x86 AVX2 SIMD
intrinsics, so the adapter only compiles on **x86_64**. On aarch64/arm64 hosts
(e.g. Apple Silicon) it can be built and run via x86_64 emulation (Rosetta
inside a Linux VM). CMake auto-detects the architecture and silently skips the
`reference_adapter` target when the reference repository is not present or the
host is non-x86_64.

```sh
# Assumes the reference is cloned next to this project at ../bamboofilters
scripts/run_reference.sh all      # setup + build + synthetic + ecoli sweeps
```

This produces `results/reference_synthetic.csv` and `results/reference_ecoli.csv`.
To diff against our results:

```sh
python3 scripts/compare.py \
    --ours      results/synthetic.csv \
    --reference results/reference_synthetic.csv
```

```sh
python3 scripts/compare.py \
    --ours      results/ecoli.csv \
    --reference results/reference_ecoli.csv
```

The script prints a per-test-case table with `× = ours / reference` ratios for
insert time and peak memory. Ratios below 1 mean our implementation is faster
or uses less memory; the project target is to stay within 2× of the reference
on every configuration.

## Repository layout

```
src/             bamboo filter core (hash, segment, bamboo_filter)
benchmark/       CLI entry point, FASTA / k-mer / generator utilities
test/            unit tests (one binary per module)
scripts/         workload, plotting and comparison scripts
data/            input data (E. coli reference genome)
results/         CSV outputs and PNG figures from benchmark runs
extern/          adapter that wraps the reference Bamboo implementation
LICENSE          project licence (MIT)
LICENSE-nthash   bundled ntHash licence (MIT, by Hamid Mohamadi)
```

## Licence

This project is released under the MIT licence (see `LICENSE`).

It bundles `nthash.hpp` v1.0.5 by Hamid Mohamadi (BC Cancer Genome Sciences
Centre) for k-mer hashing. The original MIT licence is preserved in
`LICENSE-nthash`.

## References

> [1] H. Wang, H. Dai, M. Li, J. Yu, R. Gu, J. Zheng, G. Chen. *Bamboo Filters:
> Make Resizing Smooth.* Proceedings of IEEE 38th International Conference on
> Data Engineering (ICDE), 979–991, 2022.

> [2] H. Wang et al. *Bamboo Filters: Make Resizing Smooth and Adaptive.*
> IEEE/ACM Transactions on Networking, 32(5), 3776–3791, 2024.

> [3] B. Fan, D. G. Andersen, M. Kaminsky, M. D. Mitzenmacher. *Cuckoo Filter:
> Practically Better Than Bloom.* Proceedings of ACM CoNEXT, 75–88, 2014.

> [4] H. Mohamadi, J. Chu, B. P. Vandervalk, I. Birol. *ntHash: recursive
> nucleotide hashing.* Bioinformatics, 32(22), 3492–3494, 2016. Repository:
> <https://github.com/bcgsc/ntHash>
