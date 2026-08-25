#!/usr/bin/env python3
"""Regenerate the synthetic 46,XY test BAM fixture (synth_xy.bam).

Fabricates alignments (no real sequence — CIGAR-only, SEQ/QUAL '*') tiling the
autosomes at ~4x and chrX-nonPAR / chrY-euchromatin at ~2x, i.e. a 46,XY dosage
(Rx ~= Ry ~= 0.5). Chromosome lengths are realistic (chrX 155 Mb, chrY 57 Mb) so
the PAR / Y-euchromatin coordinate masks in rx_sex apply as they do on real data.

Full coverage across each chromosome matters: rx_sex takes the MEDIAN over all
non-PAR chrX windows, so sparse coverage would median to 0. Tile the whole span.

Usage (needs samtools on PATH or via a container):
    python3 make_synth_bam.py > synth_xy.sam
    samtools sort -o synth_xy.bam synth_xy.sam && samtools index synth_xy.bam
"""
import random
import sys

random.seed(1)
RL = 50000
CHROMS = [("chr1", 10_000_000), ("chr2", 10_000_000), ("chr3", 10_000_000),
          ("chrX", 155_000_000), ("chrY", 57_000_000)]

lines = ["@HD\tVN:1.6\tSO:coordinate"]
for c, ln in CHROMS:
    lines.append(f"@SQ\tSN:{c}\tLN:{ln}")

reads = []


def tile(chrom, start, end, depth):
    n = int(depth * (end - start) / RL)
    for _ in range(n):
        pos = random.randint(start + 1, max(start + 1, end - RL))
        reads.append((chrom, pos))


for c, ln in CHROMS[:3]:
    tile(c, 0, ln, 4)                 # autosomes ~4x
tile("chrX", 0, 155_000_000, 2)       # chrX ~2x (rx_sex drops PAR windows itself)
tile("chrY", 2_650_000, 26_000_000, 2)  # chrY euchromatin ~2x

reads.sort(key=lambda r: (r[0], r[1]))
for i, (c, p) in enumerate(reads):
    lines.append(f"r{i}\t0\t{c}\t{p}\t60\t{RL}M\t*\t0\t0\t*\t*")

sys.stdout.write("\n".join(lines) + "\n")
