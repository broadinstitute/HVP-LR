#!/usr/bin/env python3
"""Compute (Rx, Ry) from mosdepth *.summary.txt* and emit a calibrate.py manifest.

NOTE: summary Ry uses WHOLE-chromosome chrY mean, which includes the Yq12 DYZ1/DYZ2
satellite mass — so it UNDERESTIMATES the euchromatin-only RY_UNIT the method wants.
Rx (primary axis) is unaffected (chrX has no comparable satellite). RY here is
PROVISIONAL; recompute from regions_bed windows for the real RY_UNIT.
"""
import glob
import os
import statistics as st

# CEU trio: ground-truth karyotypes, also self-confirmed by depth (chrY≈0 for XX).
LABELS = {
    "NA12878": ("46,XX", 2, 0, "anchor"),
    "NA12891": ("46,XY", 1, 1, "anchor"),   # non-HG002 male -> valid RY anchor
    "NA12892": ("46,XX", 2, 0, "anchor"),
}
AUTO = {f"chr{i}" for i in range(1, 23)}

def parse(path):
    means = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            if f[0].endswith("_region"):
                continue
            means[f[0]] = float(f[3])
    auto = st.median([means[c] for c in AUTO if c in means])
    rx = means["chrX"] / auto
    ry = means.get("chrY", 0.0) / auto
    return auto, rx, ry

rows = ["sample_id\tn_x\tn_y\tplatform\try_role\trx\try"]
print(f"{'sample':10s} {'auto':>6s} {'Rx':>6s} {'Ry(prov)':>9s}  label")
for path in sorted(glob.glob("controls/summary/*.summary.txt")):
    sid = os.path.basename(path).split(".")[0]
    if sid not in LABELS:
        continue
    label, nx, ny, role = LABELS[sid]
    auto, rx, ry = parse(path)
    print(f"{sid:10s} {auto:6.2f} {rx:6.3f} {ry:9.3f}  {label}")
    rows.append(f"{sid}\t{nx}\t{ny}\tshort_read\t{role}\t{rx:.4f}\t{ry:.4f}")

open("controls/manifest_ceu.tsv", "w").write("\n".join(rows) + "\n")
print("\nwrote controls/manifest_ceu.tsv")
