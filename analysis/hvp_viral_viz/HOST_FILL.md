# Host species fill — alternatives + chosen strategy

## Problem

`annotated_hits.tsv` carries no host info. BFVD reference tables
(`bfvd/`) carry virus taxonomy but no host either. Need
`taxid → host_group` mapping for UMAP-by-host plot.

## Definition of `host_group`

Coarse buckets that compress 23,614 viral taxids into ~10 categories
useful for visual stratification:

| host_group | examples |
|---|---|
| `human` | HIV, HSV, HPV, SARS-CoV-2, common-cold coronaviruses |
| `vertebrate_nonhuman` | bat coronaviruses, avian influenza, fish iridoviruses |
| `arthropod` | insect baculoviruses, mosquito-vectored viruses (vector not host of replication) |
| `plant` | TMV, plant geminiviruses |
| `fungus` | mycoviruses |
| `protist` | algal viruses, mimivirus host *Acanthamoeba* |
| `bacteria` | bacteriophages (Caudoviricetes, Microviridae, etc.) |
| `archaea` | archaeal viruses |
| `metazoan_other` | nematodes, molluscs (abalone herpesvirus) |
| `unknown` | uncharacterized environmental virus |

`host_organism` retained as free-text alongside (e.g. "Homo sapiens",
"Gordonia bronchialis") for hover-text in interactive plots.

## Source ranking — best-first

Use a layered fill: each layer fills any taxid the previous left
unknown, with `host_source` recording which layer landed each value.

### 1. ICTV Virus Metadata Resource (VMR) — primary

- Single authoritative table per ICTV release (MSL40 current as of 2026).
- Columns include `Host source` (human-curated) at species / genus / family granularity.
- Pull URL: `https://ictv.global/vmr/current` → CSV → ~17k rows mapping ICTV taxon → host source.
- **Limit:** ICTV taxonomy ≠ NCBI taxid. Need ICTV name → NCBI taxid bridge.
  Two paths:
  - Use ICTV's `Genome composition` + `Virus name(s)` to fuzzy-match to NCBI taxonomy names (provided in BFVD lineage already).
  - Use `ictv-taxonomy-id` field in VMR (if released after MSL38, which carries NCBI taxid for ~80% of entries).
- Expected coverage: 60-75% of BFVD taxids hit ICTV exactly at species or genus level.

### 2. UniProt REST `organism_host` — secondary

- Each viral UniProt entry has an `Organism host` cross-reference (when
  characterized). Query API at
  `https://rest.uniprot.org/uniprotkb/{accession}.json?fields=organism_host`.
- Coverage on BFVD: estimated 30-50% (many BFVD entries are uncharacterized
  metagenome predictions where UniProt has no host annotation).
- For the unique 347k uniprots in BFVD, this is ~70k REST calls if
  batched (UniProt allows 500 ids per stream/search query → ~700 calls
  with rate-limit pause). Cacheable as `refs/uniprot_host.parquet`.

### 3. NCBI Taxonomy host field — tertiary

- NCBI taxonomy entries carry `host` in the dataset (efetch from
  `taxonomy` db). Coverage similar to UniProt but for the taxid (not
  uniprot) — fills viruses where the species is known but the UniProt
  entry is bare.
- ~24k taxids → 24k efetch calls (cacheable; batched in groups of 200).

### 4. Rule-based by lineage — quaternary

- Sweep remaining unknowns using BFVD lineage we already have.
  Mapping rules from clade → host_group:

  ```
  class Caudoviricetes        → bacteria
  family Microviridae         → bacteria
  family Inoviridae           → bacteria
  family Tectiviridae         → bacteria
  family Plasmaviridae        → bacteria
  family Corticoviridae       → bacteria
  family Sphaerolipoviridae   → bacteria / archaea (mostly bacteria; flag dual)
  family Fuselloviridae       → archaea
  family Rudiviridae          → archaea
  family Lipothrixviridae     → archaea
  family Globuloviridae       → archaea
  family Ampullaviridae       → archaea
  family Bicaudaviridae       → archaea
  family Salterprovirus       → archaea
  order Geplafuvirales        → plant
  family Geminiviridae        → plant
  family Caulimoviridae       → plant
  family Bromoviridae         → plant
  family Closteroviridae      → plant
  family Tymoviridae          → plant
  family Tobamoviridae        → plant
  family Virgaviridae         → plant
  family Mimiviridae          → protist (Acanthamoeba)
  family Marseilleviridae     → protist
  family Phycodnaviridae      → protist (algae)
  family Pithoviridae         → protist
  family Asfarviridae         → vertebrate_nonhuman (mostly suid)
  family Iridoviridae         → vertebrate_nonhuman + arthropod (flag mixed)
  family Baculoviridae        → arthropod
  family Polydnaviridae       → arthropod
  family Nudiviridae          → arthropod
  family Hytrosaviridae       → arthropod
  family Hytrosaviridae       → arthropod
  family Partitiviridae       → fungus / plant (flag mixed)
  family Totiviridae          → fungus / protist (flag mixed)
  family Hypoviridae          → fungus
  family Polyomaviridae       → vertebrate_nonhuman (mostly mammals; human via taxid override)
  family Papillomaviridae     → vertebrate_nonhuman (mostly mammals; human via taxid override)
  family Herpesviridae        → vertebrate_nonhuman (human via taxid override)
  family Coronaviridae        → vertebrate_nonhuman (human via taxid override)
  family Flaviviridae         → vertebrate_nonhuman + arthropod (vector; flag)
  family Filoviridae          → vertebrate_nonhuman (human via taxid override)
  family Retroviridae         → vertebrate_nonhuman (human via taxid override)
  ```

- After-rule sweep, force any taxid in NCBI `Viruses_human` curated
  list → `host_group = human` (overrides family default).

### 5. Scientific-name pattern match — last resort

- BFVD's scientific name often embeds the host: `Gordonia phage Dardanus`,
  `Abalone herpesvirus`, `Acanthamoeba polyphaga mimivirus`. Regex:
  - `phage` or `prophage` → bacteria
  - First word is a known prokaryote genus → bacteria (via NCBI taxonomy
    dump cross-check)
  - First word ends in `virus` and is preceded by a host species name →
    derive host_organism from preceding tokens
- This catches an estimated 20-30% of remaining unknowns. Tagged
  `host_source = name_regex` so we can downweight in any quantitative
  analysis.

## Alternative sources NOT used (and why)

- **ViralZone (Expasy)** — hand-curated, web-only, no bulk download.
  Useful for spot-checks of high-interest viruses, not bulk fill.
- **NCBI Virus** (the portal, not the taxdump) — overlaps with NCBI
  taxonomy host field; portal API is rate-limited and intended for
  interactive use.
- **CHVD / IMG/VR / GOV2** — environmental virus catalogs with host
  prediction. Useful if our novel hits cluster around uncharacterized
  BFVD entries; would attach via target uniprot → environmental virus
  catalog ID → predicted host. Out of scope for v1.
- **CRISPR-spacer-based host prediction** (e.g. iPHoP, PHIST) — strongest
  predictor for phages whose host has CRISPR records. Worth running on
  the BFVD subset that fails layers 1-5; deferred to a separate task
  since it requires downloading several large reference DBs.
- **Sequence-similarity-based host transfer** (BLAST our query against
  a host-tagged virus DB, transfer host of nearest neighbor) — circular
  for our use case (we'd be using the same hits we're trying to label).

## Cost of layers 1-3

| Layer | API calls | Wall time (single-thread) | Cache reuse |
|---|---|---|---|
| 1 ICTV VMR | 1 (CSV download) | 5 s | Forever |
| 2 UniProt REST | ~700 batched queries | ~30 min (with rate-limit) | Per uniprot |
| 3 NCBI efetch | ~120 batches of 200 | ~10 min | Per taxid |

All caches stored under `refs/cache/` as parquet; only refetch on
explicit `--refresh-cache`.

## Expected final coverage

| Source layer | Estimated cumulative coverage |
|---|---|
| ICTV VMR | 60-75% |
| + UniProt host | 75-85% |
| + NCBI taxonomy host | 80-88% |
| + Rule-based by lineage | 92-97% |
| + Name regex | 95-99% |
| `unknown` residual | 1-5% |

All `unknown` entries plot as gray in the host UMAP and are excluded
from host-stratified abundance bars (counted separately).

## Audit trail

`host_fill_audit.tsv` columns:
```
taxid, scientific_name, family, genus, host_group, host_organism,
host_source, host_confidence
```

`host_confidence` ∈ {`high`, `medium`, `low`}:
- `high` = layers 1-3 (curated)
- `medium` = layer 4 (rule by lineage)
- `low` = layer 5 (name regex)

Plots in `out/` can be regenerated against `host_confidence ≥ medium`
to drop name-regex entries when reviewer wants conservative figures.

## Chosen v1 scope

Implement layers 1, 4, 5 first (zero cost, fastest signal). Add
layers 2 + 3 once layers 1+4+5 leave a residual we can quantify. This
keeps `make refs` fast for iteration and avoids burning UniProt /
NCBI quota during dev.
