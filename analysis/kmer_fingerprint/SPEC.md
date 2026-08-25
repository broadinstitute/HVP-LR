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
- **Sites: the Broad/GATK–Picard haplotype database** (well-vetted, LD-aware,
  stood the test of time). Used only to define *which* SNPs to fingerprint.
- **Tool: `ntsm`** (Justin Chu) — counts allele k-mers from raw reads, compares
  counts to compute P(same origin). Works on ILMN + HiFi. Repo:
  https://github.com/JustinChu/ntsm
- **Reference for marker construction: hg38** (offline, one-time — to pull ref
  allele + flanks per SNP). Downstream data need not be hg38.
- **k-mer params: ntsm defaults `k=19`, `w=31`, `n ≤ 12`.** NOTE: GATK/Picard
  fingerprinting is genotype-likelihood-based and has no k-mer size; k belongs to
  ntsm. 21 is only a generic default — use ntsm's tuned 19/31.

## Fit for ILMN + HiFi (host-present, this phase)

- Both platforms ~0.1–1%, substitution-dominated → per-base error is a non-issue
  for k-mers here (the concern that would sink noisy ONT/CLR does not apply).
- Host-present → ample human coverage at markers → no coverage-floor worry now.
- **Host-depletion handling is explicitly deferred** to a later phase; it affects
  only the *match* step (a matchability gate on residual human marker coverage),
  not the *compute* step — so nothing here needs rework when we add it.

## One-time setup (build the marker file)

1. Broad haplotype DB SNPs → sites VCF on **hg38**.
2. `ntsm` marker generation from that VCF + hg38 FASTA → ref/alt k-mer FASTAs.
3. **Uniqueness QC**: keep only SNPs whose ref *and* alt k-mers are genome-unique
   (drop repeat/segdup sites — non-unique k-mers cause false hits). Record how
   many survive.
4. **Site-count check** (open risk): ntsm's tuned regime is ~10⁵ sites; the Broad
   DB is ~thousands. For binary identity this is very likely sufficient (identity
   LOD saturates with far fewer sites than relatedness estimation), but confirm by
   validation (below). If power is marginal, supplement with ntsm's shipped
   ~10⁵-site human panel — but keep the Picard sites as the vetted core.
5. **Freeze + version** the marker file. Every sample must use the same versioned
   markers to be comparable; the version is part of the fingerprint's identity.

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

- Exact Broad haplotype-DB artifact + its post-uniqueness-QC site count (drives
  the site-count risk above).
- Whether to supplement with ntsm's shipped panel if identity power is marginal.
- Storage location/format for the per-sample fingerprint artifact in the pipeline.

## Deferred (not this phase)

- Host-depletion matchability gate (per-sample residual-human-coverage + per-pair
  shared-marker floor; report low-signal as *unmatchable*, never *mismatched*).
- WDL + Docker (mirror x_dosage's scatter/gather shape:
  per-sample `ntsm count` scatter → `ntsm eval` gather → assignment/report).
