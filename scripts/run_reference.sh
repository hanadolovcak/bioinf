#!/usr/bin/env bash
# Build and run the reference Bamboo Filter implementation
# (wanghanchengchn/bamboofilters) on the same workloads as ours.
#
# Expectations:
#   - The reference is cloned next to this repo as a sibling directory, i.e.
#     <some-parent>/bioinf and <some-parent>/bamboofilters (the default for
#     `scripts/run_reference.sh setup` below).
#   - The target machine is x86_64 with AVX2 — the reference uses AVX2
#     intrinsics in src/bamboofilter/segment.hpp and will not compile without.
#
# Usage:
#   scripts/run_reference.sh setup            # clone the reference repo
#   scripts/run_reference.sh build            # configure + build (this project's
#                                              CMake auto-detects and adds the
#                                              reference_adapter target)
#   scripts/run_reference.sh synthetic        # run the matching sweep
#   scripts/run_reference.sh ecoli            # E. coli sweep
#   scripts/run_reference.sh all              # setup → build → synthetic → ecoli

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PARENT_DIR="$(cd "$ROOT_DIR/.." && pwd)"
REF_DIR_DEFAULT="$PARENT_DIR/bamboofilters"
REF_DIR="${REFERENCE_DIR:-$REF_DIR_DEFAULT}"
BIN="$ROOT_DIR/build/reference_adapter"
RESULTS_SYNTH="$ROOT_DIR/results/reference_synthetic.csv"
RESULTS_ECOLI="$ROOT_DIR/results/reference_ecoli.csv"

cmd="${1:-help}"

setup() {
    if [ -d "$REF_DIR/.git" ]; then
        echo "Reference repo already present at $REF_DIR"
    else
        git clone --depth 1 https://github.com/wanghanchengchn/bamboofilters.git "$REF_DIR"
    fi
}

build() {
    cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build" -DCMAKE_BUILD_TYPE=Release \
        -DREFERENCE_DIR="$REF_DIR"
    cmake --build "$ROOT_DIR/build" -j 4
    if [ ! -x "$BIN" ]; then
        echo
        echo "reference_adapter binary missing. The likely cause is:"
        echo " - You are on $(uname -m); the reference impl requires x86_64."
        echo " - The reference repo at $REF_DIR is missing — run 'setup' first."
        exit 1
    fi
}

run_synthetic() {
    mkdir -p "$ROOT_DIR/results"
    rm -f "$RESULTS_SYNTH"
    local data_sizes=(1000 10000 100000 1000000 10000000)
    local ks=(10 20 50 100 200)
    for n in "${data_sizes[@]}"; do
        for k in "${ks[@]}"; do
            if [ "$k" -ge "$n" ]; then continue; fi
            echo ">>> reference synthetic n=$n k=$k"
            "$BIN" --mode synthetic --data-size "$n" --k "$k" \
                --fpr-queries 10000 --seed 42 --output "$RESULTS_SYNTH"
        done
    done
    echo "wrote $RESULTS_SYNTH"
}

run_ecoli() {
    local fasta="$ROOT_DIR/data/ecoli.fa"
    if [ ! -s "$fasta" ]; then
        echo "Missing $fasta — commit it to data/ first." >&2
        exit 1
    fi
    mkdir -p "$ROOT_DIR/results"
    rm -f "$RESULTS_ECOLI"
    for k in 10 20 50 100 200; do
        echo ">>> reference ecoli k=$k"
        "$BIN" --mode ecoli --fasta "$fasta" --k "$k" \
            --fpr-queries 10000 --seed 42 --output "$RESULTS_ECOLI"
    done
    echo "wrote $RESULTS_ECOLI"
}

case "$cmd" in
    setup)     setup ;;
    build)     build ;;
    synthetic) run_synthetic ;;
    ecoli)     run_ecoli ;;
    all)       setup; build; run_synthetic; run_ecoli ;;
    help|*)
        echo "Usage: $0 {setup|build|synthetic|ecoli|all}"
        echo "Reference repo expected at \$REFERENCE_DIR or $REF_DIR_DEFAULT"
        ;;
esac
