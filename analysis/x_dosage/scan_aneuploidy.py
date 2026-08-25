#!/usr/bin/env python3
"""Find sex-chromosome aneuploidy candidates by dosage-vs-reported-sex discordance.

For every mosdepth *.summary.txt, compute whole-chromosome dosage:
    Rx = mean(chrX) / median(chr1..22)      Ry = mean(chrY) / median(chr1..22)
Join reported sex (igsr_samples.tsv). Euploid expectation:
    male   : Rx ~ 0.5,  Ry ~ 0.3-0.4 (whole-chrom, Yq12-diluted)
    female : Rx ~ 1.0,  Ry ~ 0
A sample whose dosage contradicts its reported sex is an aneuploidy candidate.

NOTE: Rx/Ry here are WHOLE-chromosome (summary-based) — deliberately coarse, only
to SHORTLIST. Confirm each hit with rx_sex.py on its regions_bed (euchromatin Ry,
bootstrap). This finds candidates from real depth; it does not fabricate a list.

Thresholds are dosage midpoints between euploid and the nearest aneuploid class,
not tuned to any expected count.
"""
import glob
import os
import statistics as st
import sys

AUTO = {f"chr{i}" for i in range(1, 23)}


def dosage(path):
    m = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            if f[0].endswith("_region"):
                continue
            m[f[0]] = float(f[3])
    auto = st.median([m[c] for c in AUTO if c in m])
    if auto <= 0:
        return None
    return m.get("chrX", 0.0) / auto, m.get("chrY", 0.0) / auto


def classify_candidate(sex, rx, ry):
    """Return (implied_karyotype, reason) if discordant with reported sex, else None."""
    if sex == "male":
        # euploid male ~ (0.5, 0.3-0.4)
        if rx > 1.30:
            return "48,XXXY?", f"male Rx={rx:.2f} (~3X)"
        if rx > 0.75:
            return "47,XXY", f"male Rx={rx:.2f} (~2X)"
        if ry > 0.60:
            return "47,XYY", f"male Ry={ry:.2f} (~2Y)"
        if ry < 0.10:
            return "45,X / LOY", f"male Ry={ry:.2f} (~0Y)"
    elif sex == "female":
        # euploid female ~ (1.0, ~0)
        if rx > 1.30:
            return "47,XXX", f"female Rx={rx:.2f} (~3X)"
        if rx < 0.75:
            return "45,X", f"female Rx={rx:.2f} (~1X)"
        if ry > 0.15:
            return "XXY? (reported F)", f"female Ry={ry:.2f} (Y present)"
    return None


def main(summary_dir, igsr_tsv):
    sex = {}
    with open(igsr_tsv) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            sex[f[0]] = f[1].strip().lower()

    hits = []
    n = 0
    for path in sorted(glob.glob(os.path.join(summary_dir, "*.summary.txt"))):
        sid = os.path.basename(path).split(".")[0]
        d = dosage(path)
        if d is None:
            continue
        n += 1
        rx, ry = d
        s = sex.get(sid, "?")
        cand = classify_candidate(s, rx, ry)
        if cand:
            hits.append((sid, s, rx, ry, cand[0], cand[1]))

    print(f"scanned {n} samples; {len(hits)} aneuploidy candidates\n")
    print(f"{'sample':10s} {'sex':7s} {'Rx':>5s} {'Ry':>5s}  {'implied':16s} reason")
    for sid, s, rx, ry, kar, reason in sorted(hits, key=lambda h: h[4]):
        print(f"{sid:10s} {s:7s} {rx:5.2f} {ry:5.2f}  {kar:16s} {reason}")
    return hits


if __name__ == "__main__":
    sd = sys.argv[1] if len(sys.argv) > 1 else "controls/summary_all"
    ig = sys.argv[2] if len(sys.argv) > 2 else "controls/igsr_samples.tsv"
    main(sd, ig)
