# Bamboo Filter (C++17)

A from-scratch C++17 implementation of the **Bamboo Filter** approximate-membership
data structure (Wang et al., IEEE/ACM Transactions on Networking, 2024), tested on
synthetic DNA and on the *Escherichia coli* K-12 genome with k-mers of varying
length.

Project for the **Bioinformatics** course (FER, 2025/2026), task (2) — MDL track.
Link to course page:
[https://www.fer.unizg.hr/predmet/bio](https://www.fer.unizg.hr/predmet/bio).

## Authors

- Ana Angelov
- Hana Dolovčak

## Licence

This project is released under the MIT licence (see `LICENSE`).

It bundles `nthash.hpp` v1.0.5 by Hamid Mohamadi (BC Cancer Genome Sciences
Centre) for k-mer hashing. The original MIT licence is preserved in
`LICENSE-nthash`. Citation:

> Mohamadi, H., Chu, J., Vandervalk, B. P., & Birol, I. (2016).
> *ntHash: recursive nucleotide hashing.*
> Bioinformatics, 32(22), 3492–3494.
> Repository: <https://github.com/bcgsc/ntHash>

## Requirements

- Linux (Ubuntu 22.04 or similar)
- CMake ≥ 3.16
- g++ ≥ 9 *or* clang++ ≥ 10 with C++17 support
- Python 3 with matplotlib (only for plotting)

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
```

This produces:

- `build/libbamboo.a` — static library with all of `src/*`
- `build/benchmark` — CLI used to run the workloads
- `build/test_segment`, `build/test_bamboo`, `build/test_kmer`, `build/test_fasta`
  — unit tests (also runnable as a suite via `ctest`)

For a debug build with AddressSanitizer:

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j 4
```

## Running the unit tests

```sh
cd build && ctest --output-on-failure
```

## Running the benchmark

The `benchmark` executable accepts CLI flags only (no IDE, no config files):

```sh
./build/benchmark --mode synthetic \
    --data-size 1000000 --k 20 \
    --fp-bits 12 --initial-segments 8192 \
    --fpr-queries 50000 \
    --output results/synthetic.csv
```

Output is appended to the CSV named by `--output`. The header row is written
automatically the first time the file is created.

### Synthetic-data sweep (10³ … 10⁷ bases, k ∈ {10, 20, 50, 100, 200})

```sh
scripts/run_synthetic.sh
```

### E. coli genome sweep

The reference genome (NCBI accession **U00096.3**, *E. coli* K-12 substrain
MG1655, ~4.6 Mb) is committed under `data/ecoli.fa`. To run the sweep:

```sh
scripts/run_ecoli.sh
```

### Generating figures

```sh
python3 scripts/plot_results.py --csv results/synthetic.csv --prefix synthetic
python3 scripts/plot_results.py --csv results/ecoli.csv     --prefix ecoli
```

Figures are written to `results/figures/`.

## CLI reference

```
--mode {synthetic|ecoli}      workload type (default synthetic)
--k N                          k-mer length (default 20)
--data-size N                  synthetic sequence length (default 100000)
--fasta PATH                   FASTA path (required for ecoli mode)
--seed S                       RNG seed (default 42)
--fp-bits N                    fingerprint bits 1..15 (default 8)
--initial-segments N           initial bamboo segment count, must be a power of two
--bucket-per-segment N         buckets per segment, must be a power of two
--fpr-queries N                number of negative queries used to estimate FPR
--output PATH                  CSV file to append the row to
```

## Comparing with the reference implementation

The reference Bamboo Filter from
[wanghanchengchn/bamboofilters](https://github.com/wanghanchengchn/bamboofilters)
ships as a C++11 library; it does not ship a CLI that matches ours. We bundle a
small adapter (`extern/reference_adapter.cpp`) that uses the same CLI flags and
the same CSV schema as `benchmark`, but wraps *their* `BambooFilter` class.

**Platform constraint.** The reference uses x86 AVX2 intrinsics
(`src/bamboofilter/segment.hpp`), so the adapter only builds on x86_64. CMake
detects this automatically; on aarch64/arm64 (e.g. Apple Silicon) the
`reference_adapter` target is silently skipped. Run the comparison on the FER
Linux x86_64 machines.

```sh
# Assumes you cloned the reference next to bioinf, i.e. <parent>/bamboofilters
scripts/run_reference.sh all          # setup + build + synthetic + ecoli sweeps
```

This produces `results/reference_synthetic.csv` and `results/reference_ecoli.csv`
in the exact schema as our own CSVs. Then diff:

```sh
python3 scripts/compare.py \
    --ours      results/synthetic.csv \
    --reference results/reference_synthetic.csv
```

prints a row-per-test-case table with `× = ours / reference` ratios. Rows whose
ratio exceeds 2 are flagged. The FER course rubric awards points for ratio ≤ 2,
penalises –10 for ratio > 2, and –15 for ratio > 3.

## Repository layout

```
src/             — bamboo filter core (hash, segment, bamboo_filter)
benchmark/       — CLI entry point, FASTA / k-mer / generator utilities
test/            — unit tests (one binary per module)
scripts/         — workload, plotting and comparison scripts
data/            — input data (E. coli reference genome)
results/         — CSV outputs and PNG figures from benchmark runs
extern/          — adapter that wraps the reference Bamboo implementation
LICENSE          — project licence (MIT)
LICENSE-nthash   — bundled ntHash licence (MIT, by Hamid Mohamadi)
```

## References

> Wang H., Dai H., Chen S., Li M., Gu R., Chai H., Zheng J., Chen Z., Li S.,
> Deng X., Chen G. (2024). *Bamboo Filters: Make Resizing Smooth and Adaptive.*
> IEEE/ACM Transactions on Networking, 32(5), 3776–3791.

> Fan B., Andersen D. G., Kaminsky M., Mitzenmacher M. D. (2014).
> *Cuckoo Filter: Practically Better Than Bloom.*
> Proceedings of CoNEXT '14, 75–88.

> Mohamadi H., Chu J., Vandervalk B. P., Birol I. (2016).
> *ntHash: recursive nucleotide hashing.*
> Bioinformatics, 32(22), 3492–3494.
