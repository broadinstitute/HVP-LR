#!/usr/bin/env python3
"""Select a balanced euploid Tier-2 set for SIGMA_*_MODEL calibration.

For each mosdepth summary compute whole-chrom Rx, Ry; join reported sex + superpop
(igsr). Keep only clean euploids (dosage consistent with reported sex), then pick N
per (superpopulation x sex) closest to the euploid dosage center — a spread across
the 5 superpopulations that stabilizes the within-group SD (SIGMA), per
calibration_samples.md Tier 2.

Emits a gsutil fetch list of regions_bed for the picks NOT already downloaded, plus
a labeled manifest for calibrate.py.

Unrelatedness: 1kGP has trios; strict founder-filtering needs the founder list. This
selector does not enforce it (a stray relative negligibly affects a ~50-sample SD).
Pass a founder file via --founders to enforce if desired.
"""
import argparse
import glob
import os
import statistics as st

AUTO = {f"chr{i}" for i in range(1, 23)}
BUCKET = ("gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/1KGP/"
          "alignment/statistics/mosdepth/regions_bed")


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


def is_euploid(sex, rx, ry):
    if sex == "male":
        return 0.40 <= rx <= 0.65 and 0.20 <= ry <= 0.50
    if sex == "female":
        return 0.85 <= rx <= 1.15 and ry < 0.05
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-dir", default="controls/summary_all")
    ap.add_argument("--igsr", default="controls/igsr_samples.tsv")
    ap.add_argument("--per-group", type=int, default=5, help="picks per superpop x sex")
    ap.add_argument("--have-dir", default="controls/regions_1mb",
                    help="already-downloaded regions (skip re-fetch)")
    ap.add_argument("--founders", help="optional founder list to enforce unrelatedness")
    args = ap.parse_args()

    meta = {}
    with open(args.igsr) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            meta[f[0]] = (f[1].strip().lower(), f[5].strip())  # sex, superpop code

    founders = None
    if args.founders and os.path.exists(args.founders):
        founders = set(open(args.founders).read().split())

    have = {os.path.basename(p).split(".")[0]
            for p in glob.glob(os.path.join(args.have_dir, "*.regions.bed.gz"))}

    # group -> list of (dist_to_center, sample, sex, superpop, rx, ry)
    groups = {}
    for path in sorted(glob.glob(os.path.join(args.summary_dir, "*.summary.txt"))):
        sid = os.path.basename(path).split(".")[0]
        if sid not in meta:
            continue
        sex, sp = meta[sid]
        if founders is not None and sid not in founders:
            continue
        d = dosage(path)
        if d is None:
            continue
        rx, ry = d
        if not is_euploid(sex, rx, ry):
            continue
        center = 1.0 if sex == "female" else 0.5
        groups.setdefault((sp, sex), []).append((abs(rx - center), sid, sex, sp, rx, ry))

    picks = []
    for key in sorted(groups):
        for _, sid, sex, sp, rx, ry in sorted(groups[key])[: args.per_group]:
            picks.append((sid, sex, sp, rx, ry))

    # manifest
    with open("controls/manifest_tier2.tsv", "w") as mf:
        mf.write("sample_id\tn_x\tn_y\tplatform\try_role\tregions\n")
        for sid, sex, _sp, _rx, _ry in picks:
            nx, ny = (2, 0) if sex == "female" else (1, 1)
            mf.write(f"{sid}\t{nx}\t{ny}\tshort_read\tanchor\t"
                     f"controls/regions_1mb/{sid}.regions.bed.gz\n")

    new = [p for p in picks if p[0] not in have]
    print(f"selected {len(picks)} euploids "
          f"({sum(1 for p in picks if p[1]=='female')}F / "
          f"{sum(1 for p in picks if p[1]=='male')}M) across "
          f"{len({p[2] for p in picks})} superpops; {len(new)} need fetching\n")
    for key in sorted(groups):
        n = sum(1 for p in picks if (p[2], p[1]) == key)
        print(f"  {key[0]:4s} {key[1]:6s}: {n}")

    ids = " ".join(sorted(p[0] for p in new))
    print("\n# fetch (requester-pays; set P=YOUR_PROJECT):")
    print(f"gsutil -u $P -m cp {BUCKET}/{{{','.join(sorted(p[0] for p in new))}}}"
          ".regions.bed.gz controls/regions_bed/")
    open("controls/tier2_fetch_ids.txt", "w").write(ids + "\n")


if __name__ == "__main__":
    main()
