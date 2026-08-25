# x_dosage — sex-chromosome karyotyping for sample-swap QC

Coverage-based sex-chromosome **karyotype** caller (46,XX / 46,XY / 47,XXY /
47,XYY / 47,XXX / 45,X / 48,XXYY / **OTHER**) on **T2T-CHM13v2.0**. Part of the
HVP sample-swap fingerprinting effort.

> **⚠️ Maintenance moved.** The scripts + calibrated configs are now maintained in
> **`docker/x_dosage/src/`** and baked into the **`x_dosage` image**
> (`…/hvp-longread-containers/x_dosage:<VERSION>`). Run them via the WDL
> **`wdl/pipelines/TechAgnostic/QC/SexChromKaryotype.wdl`** (mosdepth → classify),
> not from here. This `analysis/` directory is the **frozen calibration &
> validation record** (design docs + the calibration/aneuploidy provenance);
> the code that ships lives in the image. Edit scripts there, bump the image
> `VERSION`, and update the WDL `tool_version` to match.

## Why

HVP is a viral-genomics study whose samples carry **human host reads** (not
host-depleted → normal human coverage), in **both short- and long-read** data,
with **no SNP array**. To catch sample swaps / mislabels we fingerprint each
sample two ways:

1. **Sex-chromosome karyotype** (this folder) — a cheap, reference-light per-
   sample call from coverage alone. A karyotype that disagrees with the manifest
   is a swap flag.
2. Identity fingerprint (`ntsm` / `somalier`) — separate component, not here yet.

Coverage-based karyotyping needs no variant calling and works identically on SR
and LR (only the masked coordinates are reference-specific).

## How it works

Two dosage ratios per sample, from mosdepth windows:

```
Rx = median(chrX non-PAR windows)      / median(autosome windows)   ~ n_X copies * RX_UNIT
Ry = median(chrY euchromatin windows)  / median(autosome windows)   ~ n_Y copies * RY_UNIT
```

Each candidate karyotype `(n_X, n_Y)` predicts a `(Rx, Ry)`. A Gaussian posterior
over the karyotype set **plus a diffuse OTHER class** (mosaics / contamination /
structural) gives the call and a confidence. Per-sample spread combines the
calibrated model SIGMA with a bootstrap SE over windows.

Two details that matter:
- **chrX excludes PAR1/PAR2** (they carry Y homology → dilute the X signal).
- **chrY uses only euchromatin (~2.65–26 Mb)**. Whole-chromosome Y is dominated
  by the Yq12 DYZ1/DYZ2 satellite mass, which collapses Ry even in males.
  Euchromatin-only Ry is what makes XX-vs-XXY and XY-vs-XYY separable.

## Files

**Code (now in `docker/x_dosage/`):**

| File | What |
|------|------|
| `docker/x_dosage/src/rx_sex.py` | The classifier. `--config <toml>`; stdlib-only core. |
| `docker/x_dosage/src/calibrate.py` | Labeled controls → per-platform config (`*_UNIT`, `SIGMA_*`, dosage-linearity check). |
| `docker/x_dosage/src/aggregate_to_1mb.py` | Re-bin fine (500 bp) mosdepth windows up to 1 Mb (only needed for pre-binned calibration data; the WDL runs mosdepth `--by 1000000` directly). |
| `docker/x_dosage/src/scan_aneuploidy.py` | Find aneuploidy candidates by dosage-vs-reported-sex discordance over mosdepth summaries. |
| `docker/x_dosage/src/select_tier2.py` | Pick a balanced euploid panel (5 superpops × N/sex) to stabilize SIGMA. |
| `docker/x_dosage/src/config.example.toml` | Format template = the built-in defaults. |
| `docker/x_dosage/src/configs/short_read.toml` | **Calibrated** short-read config (baked into the image; the WDL default). |
| `docker/x_dosage/tests/` | pytest: synthetic clean/mosaic/LOY + calibrate unit tests. |
| `docker/x_dosage/{Dockerfile,Makefile,env.yaml,…}` | The six-file image build. |

**Pipeline (WDL):**

| File | What |
|------|------|
| `wdl/pipelines/TechAgnostic/QC/SexChromKaryotype.wdl` | Workflow: BAM → mosdepth → karyotype call + confidence. |
| `wdl/tasks/QC/Mosdepth.wdl` | Windowed depth (1 Mb, MAPQ≥20). |
| `wdl/tasks/QC/RxSexKaryotype.wdl` | Classify via `rx_sex.py`; emits typed call/confidence/Rx/Ry. |

**This directory (frozen record):**

| File | What |
|------|------|
| `HANDOFF.md`, `calibration_samples.md` | Design rationale + calibration-set spec. Read these first. |
| `pipeline.sh` | Reference glue (FASTQ → align → mosdepth → classify) the WDL was ported from. |
| `controls/` | Downloaded depth data + per-run configs (gitignored — local only). |

## Setup

```bash
uv sync            # env with ruff + pytest
uv run pytest -q   # 18 tests
uv run ruff check .
```

The core (`rx_sex.py`, `calibrate.py`) is stdlib-only (Python ≥3.11 for
`tomllib`) — it runs without the env; `uv` is only for tests/lint.

## Usage

**Classify** (mosdepth `*.regions.bed(.gz)`, T2T-CHM13v2.0, run with
`--by 1000000 -Q 20`):

```bash
python rx_sex.py --config configs/short_read.toml sample.regions.bed.gz
# sample  auto_dp  Rx  Rx_CI  Ry  Ry_CI  call  conf  runner_up
```

If your mosdepth used finer windows (e.g. 500 bp), aggregate first:

```bash
python aggregate_to_1mb.py sample.500bp.regions.bed.gz sample.1mb.regions.bed.gz
```

**Calibrate** a new platform from labeled controls (manifest columns:
`sample_id  n_x  n_y  platform  ry_role  regions` — or precomputed `rx`/`ry`):

```bash
python calibrate.py --manifest controls/manifest_euploid56.tsv --out-dir configs/
```
`ry_role=ceiling` marks a reference-Y sample (e.g. HG002) so it is excluded from
the RY_UNIT anchor (it maps to its own Y → biased high). Anchor RY on non-ref
males.

**Find aneuploidy test samples** (dosage discordant with reported sex):

```bash
python scan_aneuploidy.py controls/summary_all controls/igsr_samples.tsv
```

## Calibration status

Calibrated + validated on **real 1kGP data aligned to T2T-CHM13v2.0** (from the
AnVIL_T2T bucket; per-sample mosdepth is precomputed there — no realignment
needed).

- **Short read** (`configs/short_read.toml`): 56 euploids (5 superpopulations ×
  ~5/sex Tier-2 + CEU/YRI trios).
  `RX_UNIT=0.489  RY_UNIT=0.490  SIGMA_RX=0.0090  SIGMA_RY=0.0066`.
  X-linearity ratio 0.995 (PASS). **Euploid validation: 56/56, min conf 0.994.**
- **Aneuploidy spot-check** (real 1kGP): clean-multiple 47,XXY / 47,XYY / 45,X
  named correctly; off-multiple / mosaic cases (an XYY at Ry 0.92, a mosaic XXX
  at Rx 1.38) → OTHER by design; mosaic-X → OTHER; full loss-of-Y male → 45,X.
- **Long read**: not yet calibrated. Euploid units are LR-calibratable (GIAB);
  **aneuploidy classes stay dosage-extrapolated on the LR path** — there is no
  public long-read aneuploidy data to fit or validate. Known, preserved limitation.

## Design decisions (settled — see HANDOFF.md, do not relitigate)

- Reference: **T2T-CHM13v2.0** (`chm13v2.0_maskedY.rCRS`).
- **Rx is primary; Ry advisory.** X-dosage carries the call; Y refines it.
- **Euchromatin-only Ry**; **OTHER class is deliberate**.
- **Euploid SIGMA is intentionally tight** (~0.009 at 30×; adding samples
  tightened it). Off-multiple/mosaic aneuploidies therefore land in OTHER rather
  than being force-fit to a karyotype — for swap QC, OTHER still flags the
  anomaly. Do not widen SIGMA to name more aneuploidies without revisiting this
  (ratified 2026-08; see `rx_sex.classify` docstring).
- **Full LOY reads as 45,X** and true mosaics read as XY↔Turner intermediates —
  flag for review; don't let a low 45,X prior silently absorb them.

## Data source

Per-sample mosdepth output for 1kGP-on-CHM13v2.0 lives in the AnVIL_T2T workspace
bucket (`gs://fc-47de7dae-…/1KGP/alignment/statistics/mosdepth/`):
`summary/` (per-chrom means, tiny), `regions_bed/` (500 bp windows, ~47 MB/sample).
It is a **requester-pays** bucket — `gsutil -u <your-project> cp …`. Reported sex
+ superpopulation: `1KGP/resources/igsr_samples.tsv`; unrelated founders:
`founder_samples.txt`.

## Next

- **DONE: WDL + Docker.** `x_dosage` image + `SexChromKaryotype.wdl`
  (mosdepth → classify) verified end-to-end (synthetic 46,XY and real 1kGP XXY).
- Short-read Y-linearity rung is INSUFFICIENT (no XYY among the *euploid*
  calibration samples; confirmed separately via a Tier-3 XYY). Add one XYY to the
  euploid manifest to close it.
- Long-read euploid calibration (GIAB SR+LR matched) → emit `long_read.toml`.
- Deferred (HANDOFF TODOs): empirical Y-mask refinement (#4), SR-vs-LR
  concordance → ntsm/somalier identity layer (#6, see `../kmer_fingerprint/`).
