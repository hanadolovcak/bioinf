#!/usr/bin/env python3
"""Compare our bamboo filter results against the reference implementation."""
import argparse
import csv
import sys
from collections import defaultdict


def read(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def key(r):
    return (r["mode"], int(r["data_size"]), int(r["k"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--reference", required=True)
    args = ap.parse_args()

    ours = {key(r): r for r in read(args.ours)}
    ref = {key(r): r for r in read(args.reference)}

    print(f"{'mode':12s} {'data_size':>10s} {'k':>4s} "
          f"{'ours_ins':>10s} {'ref_ins':>10s} {'×':>5s} "
          f"{'ours_mem':>10s} {'ref_mem':>10s} {'×':>5s}")
    for k in sorted(set(ours.keys()) | set(ref.keys())):
        if k not in ours or k not in ref:
            continue
        o, r = ours[k], ref[k]
        oi, ri = float(o["insert_s"]), float(r["insert_s"])
        om, rm = float(o["peak_rss_kb"]), float(r["peak_rss_kb"])
        ti = oi / ri if ri > 0 else float("nan")
        tm = om / rm if rm > 0 else float("nan")
        warn_i = " ⚠" if ti > 2.0 else ""
        warn_m = " ⚠" if tm > 2.0 else ""
        print(f"{k[0]:12s} {k[1]:10d} {k[2]:4d} "
              f"{oi:10.4f} {ri:10.4f} {ti:5.2f}{warn_i} "
              f"{om:10.0f} {rm:10.0f} {tm:5.2f}{warn_m}")


if __name__ == "__main__":
    main()
