# k-mer fingerprint — SPEC (parked; build after x_dosage WDL)

Alignment-free per-sample genotype fingerprint for **SR↔HiFi same-individual
matching** across the overlapping ILMN and HiFi cohorts. Component 2 of the HVP
sample-provenance effort (sibling to `../x_dosage`; relates to HANDOFF TODO #6).

Status: **spec only, not built.** Settled decisions below; implement after the
x_dosage WDL/container is done.

## Goal & posture

- **Match**, don't swap-check: confirm the ILMN sample and the HiFi sample said to
  be the same individual really are (and, where labels are untrusted, assign the
  bijection). Different framing from x_dosage's within-sample swap flag.
- **Compute on ingest, use on demand.** The fingerprint is produced automatically
  when each sample is processed and stored as a first-class per-sample output. We
  may never cross-check; if we need to, the fingerprints are already there. No raw
  reads revisited later.

## Method (settled)

- **Alignment-free k-mer allele counting** at curated SNP sites. No per-sample
  alignment; the reads' own reference build is irrelevant at count time. (Decision
  B from discussion — chosen over Picard-native ExtractFingerprint, which is
  alignment-based and would fight our T2T-CHM13 world.)
- **Sites: the UNION of ntsm's shipped human panel + the Broad/GATK–Picard
  haplotype database** (see "Panel" below). ~97,180 sites.
- **Tool: `ntsm`** (Justin Chu) — counts allele k-mers from raw reads, compares
  counts to compute P(same origin). Works on ILMN + HiFi. Repo:
  https://github.com/JustinChu/ntsm — **published**: Chu et al.,
  *ntsm: an alignment-free, ultra-low-coverage, sequencing-technology-agnostic,
  intraspecies sample comparison tool for sample swap detection*, GigaScience 2024
  (doi:10.1093/gigascience/giae024). MIT, C++.
- **Reference for marker construction: hg38** (offline, one-time — to pull ref
  allele + flanks per SNP). Downstream data need not be hg38.
- **k-mer params: ntsm defaults `k=19`, `w=31`, `n ≤ 12`.** NOTE: GATK/Picard
  fingerprinting is genotype-likelihood-based and has no k-mer size; k belongs to
  ntsm. 21 (and HaplotypeCaller's assembly 10/25) are unrelated defaults — use
  ntsm's tuned 19/31. 19-mers are already genome-unique (4^19 ≫ 3×10⁹); longer k
  only adds error-fragility and neighbor-SNP collisions.

## Panel: union of ntsm + Picard sites

Measured overlap between the two candidate panels (rsID intersection, 2026-08):

| Panel | Sites | Notes |
|-------|------:|-------|
| ntsm shipped (`data/human_sites_n10.fa`) | 96,287 | prebuilt+validated k-mers; ships rotation/normalization matrices |
| Broad/Picard hg38 haplotype DB | 933 | LD-anchored, alignment-based fingerprinting panel |
| **Shared (by rsID)** | **40** | ~4% of Picard, essentially disjoint |
| **Union (planned)** | **~97,180** | 96,287 + 933 − 40 |

The panels barely overlap because they were built for different jobs: Picard's
933 are LD-anchored for *alignment-based* genotype-LOD fingerprinting; ntsm's 96k
are common-SNP-derived for *alignment-free per-SNP* k-mer counting. Picard's
"test of time" strength lives in its LD model, which ntsm's per-SNP counting does
not use — so Picard alone would be a weak, custom-built k-mer panel at ~1% of
ntsm's tested regime.

**Decision: use the UNION.** ntsm's 96k is the validated core (ready, on-regime);
the ~893 Picard-only sites are added on top so we retain compatibility/coverage
with the vetted fingerprinting SNPs at no cost to power. Build:
1. Start from ntsm's shipped 96,287-site markers (k-mers already generated + QC'd).
2. For the ~893 Picard-only rsIDs, generate ntsm-format markers with `ntsmSiteGen`
   against hg38, then apply ntsm's own uniqueness QC (bwa-align the 19-mers to
   hg38, drop any that map multiply with ≤1 mismatch; keep sites with ≥3
   non-repetitive 19-mers) — the exact procedure from the paper.
3. Merge into one versioned marker set (~97,180 sites minus any Picard adds that
   fail uniqueness QC — record the surviving count).

Caveat: ntsm's shipped **PCA rotation/normalization matrices cover only the 96k**.
They're an optional prefilter, not required for the core count + likelihood-ratio
match. Either restrict the PCA prescreen to the 96k subset, or regenerate the
matrices for the union. The identity match itself uses the full union.

Reproduce the overlap:
```
# ntsm rsIDs
curl -sL https://raw.githubusercontent.com/JustinChu/ntsm/main/data/human_sites_n10.fa \
  | grep '^>' | awk '{print $1}' | sed 's/^>//' | sort -u > ntsm_rs.txt      # 96,287
# Picard hg38 haplotype DB rsIDs
curl -sL https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt \
  | grep -oE 'rs[0-9]+' | sort -u > pic_rs.txt                                # 933
comm -12 ntsm_rs.txt pic_rs.txt | wc -l                                       # 40 shared
```

## Fit for ILMN + HiFi (host-present, this phase)

- Both platforms ~0.1–1%, substitution-dominated → per-base error is a non-issue
  for k-mers here (the concern that would sink noisy ONT/CLR does not apply).
- Host-present → ample human coverage at markers → no coverage-floor worry now.
- **Host-depletion handling is explicitly deferred** to a later phase; it affects
  only the *match* step (a matchability gate on residual human marker coverage),
  not the *compute* step — so nothing here needs rework when we add it.

## One-time setup (build the union marker file)

1. Start from ntsm's shipped 96,287-site markers (`data/human_sites_n10.fa` +
   `human_sites_center.txt`). Already generated and uniqueness-QC'd by the authors.
2. Compute the ~893 Picard-only rsIDs (Picard DB minus the 40 already in ntsm).
3. For those, `ntsmSiteGen` against **hg38** → ref/alt 19-mer markers, then ntsm's
   uniqueness QC (bwa-align 19-mers to hg38; drop any mapping multiply with ≤1
   mismatch; keep sites with ≥3 non-repetitive 19-mers). Record how many survive.
4. **Merge** ntsm 96k + surviving Picard-only markers → the union set (~97,180
   minus QC failures — record the final count).
5. **Freeze + version** the union marker file. Every sample must use the same
   versioned markers to be comparable; the version is part of the fingerprint's
   identity. (PCA prefilter matrices ship for the 96k only — see Panel caveat.)

## Per-sample task (automatic, on ingest)

- `ntsm count <markers> <sample FASTQs>` → one small per-sample fingerprint
  (allele counts per marker). Same markers/k for ILMN and HiFi.
- Emit as a stored per-sample output. Cheap (single streaming k-mer pass),
  self-contained (no other sample needed). This is the only step wired into the
  processing pipeline now.

## Cross-check (later, on demand — NOT built now)

- All-pairs `ntsm eval` over {ILMN} × {HiFi} → P(same origin) per pair.
- **Confirm mode**: each declared SR↔HiFi pair should score at the identity tier;
  flag any that don't. **De novo mode**: solve the ILMN×HiFi assignment; flag
  ambiguous (two near-tied candidates).
- Threshold at the **identity** tier, not relatedness (else 1st-degree relatives
  pass). Take the valley of the bimodal same-vs-different score distribution.
- Watch cohort relatedness; watch for the depletion matchability gate (later).

## Validation plan

- On the known SR/HiFi overlap set: same-vs-different score distributions should
  be cleanly bimodal. This both sets the threshold and proves the reduced Picard
  site count is enough. Include a known-good pair as a positive control.

## Open decisions / to confirm at build time

- **Panel: RESOLVED — union of ntsm 96k + Picard 933 (~97,180).** See Panel section.
- How many Picard-only sites survive ntsm's uniqueness QC (record the final union count).
- PCA prefilter: restrict to the 96k subset, or regenerate rotation/normalization
  matrices for the union? (Optional prescreen; not needed for the core match.)
- Storage location/format for the per-sample fingerprint artifact in the pipeline.
- Which hg38 FASTA build for `ntsmSiteGen` (must match the rsID coordinates).

## Deferred (not this phase)

- Host-depletion matchability gate (per-sample residual-human-coverage + per-pair
  shared-marker floor; report low-signal as *unmatchable*, never *mismatched*).
- WDL + Docker (mirror x_dosage's scatter/gather shape:
  per-sample `ntsm count` scatter → `ntsm eval` gather → assignment/report).
