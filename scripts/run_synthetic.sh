#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT_DIR/build/benchmark"
RESULTS="$ROOT_DIR/results/synthetic.csv"

if [ ! -x "$BIN" ]; then
    echo "Build the project first: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build" >&2
    exit 1
fi
mkdir -p "$ROOT_DIR/results"
rm -f "$RESULTS"

DATA_SIZES=(1000 10000 100000 1000000 10000000)
KS=(10 20 50 100 200)
FP_BITS=12

# pick initial_segments adaptive to the data size so the per-test fixed
# memory overhead stays in the ballpark of the reference's auto-scaled
# capacity (their constructor sizes the hash table from the expected item
# count). Without this we pay 8192 segments × ~2 kB ≈ 16 MB even for n=1000.
pick_init() {
    local n=$1
    if   [ "$n" -le 2000 ];      then echo 16
    elif [ "$n" -le 20000 ];     then echo 32
    elif [ "$n" -le 200000 ];    then echo 256
    elif [ "$n" -le 2000000 ];   then echo 2048
    else                              echo 8192
    fi
}

for n in "${DATA_SIZES[@]}"; do
    init=$(pick_init "$n")
    for k in "${KS[@]}"; do
        if [ "$k" -ge "$n" ]; then continue; fi
        echo ">>> synthetic n=$n k=$k init_segments=$init"
        "$BIN" --mode synthetic --filter bamboo \
            --data-size "$n" --k "$k" \
            --fp-bits "$FP_BITS" --initial-segments "$init" \
            --fpr-queries 10000 --seed 42 \
            --output "$RESULTS"
    done
done

echo "wrote $RESULTS"
