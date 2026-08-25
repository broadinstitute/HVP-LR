# Calibration set for the (Rx, Ry) karyotype classifier

Reference for everything below: **T2T-CHM13v2.0** (`chm13v2.0_maskedY.rCRS.fa`).
Per-sample input = mosdepth `--by 1000000 -Q 20` regions, then `rx_sex.py`.

## Tier 1 — euploid anchors (GIAB): matched short + long read, multi-ancestry

Pins RX_UNIT (1-copy from XY, 2-copy from XX), RY_UNIT, and per-platform SIGMA.
Data: Illumina + PacBio HiFi + ONT at
`ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/{AshkenazimTrio,ChineseTrio}/`
(HG001/NA12878 via IGSR / Platinum Genomes).

| Sample | Alias    | Karyotype | Ancestry   | Role   | Platforms          |
|--------|----------|-----------|------------|--------|--------------------|
| HG001  | NA12878  | 46,XX     | CEU/Utah   | pilot  | SR + HiFi + ONT    |
| HG002  | NA24385  | 46,XY     | Ashkenazi  | son    | SR + HiFi + ONT (= ref Y; mappability CEILING) |
| HG003  | NA24149  | 46,XY     | Ashkenazi  | father | SR + HiFi + ONT    |
| HG004  | NA24143  | 46,XX     | Ashkenazi  | mother | SR + HiFi + ONT    |
| HG005  | NA24631  | 46,XY     | Han Chinese| son    | SR + HiFi + ONT    |
| HG006  | NA24694  | 46,XY     | Han Chinese| father | SR + HiFi + ONT    |
| HG007  | NA24695  | 46,XX     | Han Chinese| mother | SR + HiFi + ONT    |

RY_UNIT de-bias: HG002 maps to its own Y -> ceiling. Anchor RY_UNIT on the
non-HG002 males (HG005/HG006 + 1kGP males across Y-haplogroups); treat HG002 as
upper bound.

## Tier 2 — euploid dispersion + thresholds (1kGP high-coverage, T2T-aligned)

Fully open, aligned to T2T-CHM13v2.0 via XYalign. Add ~20-30 UNRELATED euploid
per sex, spanning the 5 super-populations, to stabilize SIGMA_*_MODEL and set
decision thresholds against an empirical FPR.
Source: IGSR / 1000genomes; T2T-recalled set referenced from github.com/marbl/CHM13.

## Tier 3 — aneuploidy class validation

Read exact sample IDs from the per-sample allosome-ploidy calls in:
  Byrska-Bishop et al., Cell 2022 (3,202-sample 1kGP), supplement "Ploidy of
  allosomes" (Fig S1A / Table S1).
Filter that table for X/Y copy number, pull matching T2T CRAMs. Expect a few each:
  47,XXY  47,XYY  47,XXX  45,X    (48,XXYY rare/absent in 1kGP)
These validate class MEANS (which are dosage-predicted); no per-class SD needed.

Rare tail (48,XXYY, 49,XXXXY): EGA phs002481 (114 SCA LCLs) — access-controlled,
RNA-seq-weighted. Fallback only.

## How many (summary)

- Run at all:            2   (NA12878 XX + HG002 XY) -> both per-copy units.
- Solid SR calibration:  ~20-30 euploid / sex (1kGP) + 7 GIAB cross-platform.
- Class validation:      every XXY/XYY/XXX/45,X in the ploidy table (a few each).
- Long read:             euploid-only (7 GIAB). Aneuploidy classes stay
                         dosage-extrapolated on the LR path — a known limitation.

## Calibration procedure

1. mosdepth + rx_sex.py on each labeled sample (SR and LR pipelines SEPARATELY).
2. Per sex/platform: read the Rx and Ry cluster centers and within-group SDs
   straight off the output.
3. Set RX_UNIT, RY_UNIT = per-copy centers; SIGMA_RX_MODEL, SIGMA_RY_MODEL =
   within-group SDs. Maintain one constant set per platform.
4. Confirm dosage linearity: XX(2X) ~= 2 * XY(1X) on the X axis; XYY(2Y) ~=
   2 * XY(1Y) on the Y axis (Tier-3 XYY samples).
5. Use HG002's own chrY per-window depth profile to place the Y_EUCHROMATIN
   boundary empirically (where coverage falls into the Yq12 satellites).

## Adjacent tool
XYalign (Webster et al. 2019) — sex-chromosome-aware ploidy estimation /
remapping; used for the 1kGP T2T alignments. Check before reinventing the
coverage-normalization step.
