#!/usr/bin/env python3
"""Aggregate fine mosdepth windows (e.g. 500bp) up to 1Mb bins.

The AnVIL_T2T bucket's regions_bed uses 500bp windows (~6M/sample); rx_sex.py's
method and bootstrap are designed for --by 1000000 (calibration_samples.md). Re-bin
depth to 1Mb (length-weighted mean) so the input matches the design and the bootstrap
is tractable. Output is a plain gzipped BED with the same 4 columns mosdepth emits.
"""
import gzip
import os
import sys

BIN = 1_000_000


def flush(fh, chrom, acc):
    for b in sorted(acc):
        wsum, lsum = acc[b]
        if lsum <= 0:
            continue
        s = b * BIN
        fh.write(f"{chrom}\t{s}\t{s + BIN}\t{wsum / lsum:.4f}\n")


def rebin(src, dst):
    cur = None
    acc = {}
    with gzip.open(src, "rt") as ih, gzip.open(dst, "wt") as oh:
        for line in ih:
            f = line.split("\t")
            c, s, e, d = f[0], int(f[1]), int(f[2]), float(f[3])
            if c != cur:
                if cur is not None:
                    flush(oh, cur, acc)
                cur, acc = c, {}
            ln = e - s
            b = s // BIN
            a = acc.get(b)
            if a is None:
                acc[b] = [d * ln, ln]
            else:
                a[0] += d * ln
                a[1] += ln
        if cur is not None:
            flush(oh, cur, acc)


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    rebin(src, dst)
    print(f"rebinned {os.path.basename(src)} -> {os.path.basename(dst)}")
