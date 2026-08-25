# HANDOFF — sex-chromosome karyotyping for sample-swap QC

## What this is
Part of a **sample-swap fingerprinting** effort for a viral-genomics study with
**human host reads** (NOT host-depleted → standard host coverage), **both short and
long read**. No SNP array available. Two layers:

1. **Sex / sex-chromosome karyotype from host coverage**  ← THIS COMPONENT (working).
   First-pass swap filter + aneuploidy flag. Catches ~half of random relabels on
   its own and gives a per-sample karyotype call with a confidence.
2. **Identity fingerprint** (future): `ntsm` (alignment-free, low-cov, tech-agnostic
   variant k-mers — fits the short+long mix) as primary; `somalier` if aligning.
   Build the ntsm site panel on the SAME reference used here to avoid spurious
   mismatches. Not started.

Cross-platform tie-in: run the short-read and long-read data for a sample
separately and require the sex/karyotype calls to agree — that is a free
short-vs-long swap check before any k-mer work.

## Method (this component)
Map host reads to **T2T-CHM13v2.0**, window depth with mosdepth, classify from the
joint (Rx, Ry) signal.

    Rx = median(chrX non-PAR windows)     / median(autosome windows)  ~ n_X * RX_UNIT
    Ry = median(chrY euchromatin windows) / median(autosome windows)  ~ n_Y * RY_UNIT

Each karyotype = copy-number pair (n_X, n_Y) → predicted (Rx, Ry). `rx_sex.py`
scores the observed pair under every class + a diffuse OTHER (mosaic/contamination/
structural), returns the posterior. Per-sample sigma = sqrt(SIGMA_MODEL^2 +
bootstrap_SE^2) per axis, so low coverage widens and loses confidence automatically.

## Design decisions (do NOT relitigate)
- **Reference: `chm13v2.0_maskedY.rCRS.fa`** (UCSC-named analysis set, Y-PAR masked).
  Not the RefSeq GCF_* fasta — chrom names must be chr-prefixed for the script.
- **X-dosage (Rx) is primary; Ry is secondary/advisory.** Rx is a clean 2-fold,
  robust to low coverage, and immune to loss-of-Y. Y is repeat-/mismap-/LOY-confounded.
- **Ry uses the proximal Y EUCHROMATIN only**, not the whole non-PAR Y. Yq12 is ~34 Mb
  of DYZ1/DYZ2 satellite (>half of chrY); a median over all non-PAR Y windows collapses
  to ~0 even in males. That is why Y_EUCHROMATIN is windowed to ~2.65–26 Mb.
- **OTHER class is deliberate** — prevents overconfident forcing of off-grid samples
  (mosaics, contamination, LOY) onto the nearest karyotype. Fractional Rx/Ry → OTHER.
- **Everything is dosage-based** → calibration is the load-bearing step, not the model.

## Verified constants (T2T-CHM13v2.0)
- chrX length 154,259,566; PAR1 (0, 2,394,410); PAR2 (153,925,834, 154,259,566)
  — from official `chm13v2.0_PAR.bed`.
- chrY length 62,460,029; non-PAR Y (2,458,320, 62,122,809); Yq12 satellite ≈ 34 Mb.
- Y_EUCHROMATIN (2,650,000, 26,000,000) is APPROXIMATE round numbers — see TODO #4.
- Classes: 46,XX / 46,XY / 47,XXY / 47,XYY / 47,XXX / 45,X / 48,XXYY + OTHER.
  Extend via the KARYOTYPES table (48,XXXY=(3,1), 48,XXXX=(4,0), 49,XXXXY=(4,1)…).

## Known limitations (keep visible; don't paper over)
- Balanced/structural rearrangements that preserve dosage are invisible.
- True mosaics land in OTHER by design (fractional Rx/Ry is their signature).
- LOY in older blood → real males read as XY↔Turner intermediate. Expected; flag for
  review, don't let a low 45,X prior silently absorb them.
- **Long-read aneuploidy classes are dosage-EXTRAPOLATED** — no public long-read
  aneuploidy data exists. Only euploid units are LR-calibratable (GIAB). State it.

## TODO (where to start)
1. **Pipeline wrapper**: fastq → aligned BAM (SR: bwa-mem2; LR: minimap2) → mosdepth
   → rx_sex.py, per sample, parallel across a cohort. SR and LR pipelines separate.
   See `pipeline.sh` for the exact commands. Prefer a real workflow engine
   (Snakemake/Nextflow) over the shell script for the cohort run. Modern tooling:
   uv-managed env, ruff.
2. **Config-driven calibration**: refactor rx_sex.py to load RX_UNIT / RY_UNIT /
   SIGMA_RX_MODEL / SIGMA_RY_MODEL / class table from a per-platform config instead
   of module-level constants. Target format in `config.example.toml`. Then a small
   `calibrate.py` that ingests labeled-control outputs and emits that config
   (cluster centers → *_UNIT, within-group SD → SIGMA_*_MODEL), one config per
   platform (SR vs LR).
3. **Acquire calibration data + calibrate** — follow `calibration_samples.md`.
   Minimum: NA12878 (XX) + HG002 (XY). RY_UNIT de-bias: HG002 IS the reference Y →
   ceiling; anchor on non-HG002 males. Validate dosage linearity (2X≈2·1X, 2Y≈2·1Y).
4. **Refine Y mask**: replace the round-number Y_EUCHROMATIN with an exact interval —
   intersect against the T2T-CHM13v2.0 censat annotation BED (drop DYZ3 centromere +
   Yq12), or derive the boundary empirically from HG002's own chrY per-window depth.
5. **Formalize tests**: the synthetic-data checks used during development
   (clean XX/XY/45,X/XXY/XYY/XXX/XXYY + mosaic + LOY) should become a pytest suite.
6. **Integrate**: per-sample SR-vs-LR karyotype concordance check → then ntsm/somalier
   identity layer (component 2).

## Env / tooling
bwa-mem2, minimap2, samtools, mosdepth; python3 (rx_sex.py is stdlib-only).
rx_sex.py sets random.seed(0) for reproducible bootstrap — drop it to inspect
run-to-run CI stability. Reference-invariant math; only the masked coordinates
are T2T-specific.

## Files in this bundle
- rx_sex.py               — the classifier (current, working)
- pipeline.sh             — align → mosdepth → classify glue + optional QC one-liners
- calibration_samples.md  — sample manifest + acquisition + calibration procedure
- config.example.toml     — TARGET format for the config refactor (TODO #2); not yet wired
- HANDOFF.md              — this file
