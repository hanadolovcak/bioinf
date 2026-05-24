#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT_DIR/build/benchmark"
FASTA="$ROOT_DIR/data/ecoli.fa"
RESULTS="$ROOT_DIR/results/ecoli.csv"

if [ ! -x "$BIN" ]; then
    echo "Build the project first" >&2
    exit 1
fi
if [ ! -s "$FASTA" ]; then
    echo "Missing $FASTA — the E. coli reference genome should be committed to data/." >&2
    exit 1
fi
mkdir -p "$ROOT_DIR/results"
rm -f "$RESULTS"

for k in 10 20 50 100 200; do
    echo ">>> ecoli k=$k"
    "$BIN" --mode ecoli --filter bamboo \
        --fasta "$FASTA" --k "$k" \
        --fp-bits 12 --initial-segments 8192 \
        --fpr-queries 10000 --seed 42 \
        --output "$RESULTS"
done

echo "wrote $RESULTS"
