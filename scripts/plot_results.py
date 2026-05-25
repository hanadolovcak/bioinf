#!/usr/bin/env python3
"""Plot benchmark results from CSV files into PNG figures."""
import argparse
import csv
import os
import sys
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not available; install with: pip install matplotlib", file=sys.stderr)
    sys.exit(1)


def read_csv(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            for key in ("data_size", "k", "num_kmers", "num_unique",
                        "peak_rss_kb", "internal_bytes", "segments",
                        "positive_hits"):
                if r.get(key) is not None and r.get(key) != "":
                    r[key] = int(r[key])
            for key in ("insert_s", "lookup_s", "fpr", "load_factor"):
                if r.get(key) is not None:
                    r[key] = float(r[key])
            rows.append(r)
    return rows


def plot_metric(rows, x_key, y_key, title, xlabel, ylabel, out_path, log_x=True, log_y=True):
    by_k = defaultdict(list)
    for r in rows:
        by_k[r["k"]].append((r[x_key], r[y_key]))
    plt.figure(figsize=(7, 5))
    for k in sorted(by_k.keys()):
        pts = sorted(by_k[k])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        plt.plot(xs, ys, marker="o", label=f"k={k}")
    if log_x:
        plt.xscale("log")
    if log_y:
        plt.yscale("log")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="CSV file produced by benchmark")
    ap.add_argument("--out-dir", default="results/figures")
    ap.add_argument("--prefix", default="synthetic")
    args = ap.parse_args()

    rows = read_csv(args.csv)
    if not rows:
        print(f"no rows in {args.csv}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(args.out_dir, exist_ok=True)

    for r in rows:
        n = r.get("num_unique") or r["num_kmers"]
        r["insert_mops"] = n / r["insert_s"] / 1e6 if r["insert_s"] > 0 else 0
        r["lookup_mops"] = n / r["lookup_s"] / 1e6 if r["lookup_s"] > 0 else 0
        r["memory_mb"] = r["peak_rss_kb"] / 1024.0

    plot_metric(rows, "num_kmers", "insert_mops",
                "Insert throughput vs # k-mers",
                "# k-mers (log)", "Insert throughput (M k-mers/s)",
                os.path.join(args.out_dir, f"{args.prefix}_insert_throughput.png"),
                log_y=False)
    plot_metric(rows, "num_kmers", "lookup_mops",
                "Lookup throughput vs # k-mers",
                "# k-mers (log)", "Lookup throughput (M k-mers/s)",
                os.path.join(args.out_dir, f"{args.prefix}_lookup_throughput.png"),
                log_y=False)
    plot_metric(rows, "num_kmers", "memory_mb",
                "Peak memory vs # k-mers",
                "# k-mers (log)", "Peak RSS (MB)",
                os.path.join(args.out_dir, f"{args.prefix}_memory.png"))
    plot_metric(rows, "num_kmers", "fpr",
                "False positive rate vs # k-mers",
                "# k-mers (log)", "FPR",
                os.path.join(args.out_dir, f"{args.prefix}_fpr.png"),
                log_y=False)


if __name__ == "__main__":
    main()
