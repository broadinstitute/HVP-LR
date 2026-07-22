# hvp_viral_viz — research journal

Chronological design log for this package. Entries were extracted from the
parent HVP-LR research journal (kept outside the repo) and copied here so the
reasoning trail travels with the code. Format mirrors the parent journal:
**Context / Decision / Why / Evidence / Outcome / Follow-up**, append-only,
UTC timestamps.

If you change your mind about a previous entry, add a new entry that
references the old one ("see 2026-06-23T22:18Z — that hypothesis was
wrong because …"). Do not edit history.

---

## 2026-06-23T22:00:00Z — viral-viz: Phase 1+2 complete (HVP-0006.1_34P)

**Context.** New repo `/workspace/viral-viz` scaffolded for single-cell-style
viz of HvpViralProteinAnnotation outputs. Goal: count matrix per virus, UMAP
colored by taxonomy + host, abundance plot, plus extra tier-2 plots. User
authorized Option C (cohort scaffold + run on one sample). Source sample:
Terra submission 2da383d1-fafb-46f0-bac4-3f4f782304d0, `HVP-0006.1_34P`.

**Decision / action.** Built end-to-end pipeline:
- `build_bfvd_refs.py` → 3 parquets from BFVD tarballs
  (uniprot_taxid, taxid_lineage, uniprot_qc). 24,434 uniprots / 23,614 taxids.
- `build_host_table.py` → layered host fill: ICTV VMR MSL41
  species/genus/family + ~60 family LINEAGE_RULES + name regex + ~80 curated
  HUMAN_TAXIDS. Output: `refs/taxid_host.parquet` + audit TSV.
- `ingest.py` → m8 → AnnData on disk. DISCOVERY threshold
  `evalue≤1e-5 AND bits≥50 AND alnlen≥50`; HIGH_CONF `1e-10/300/80`.
- `plots.py` (tier 1) + `plots_tier2.py` (tier 2) → 10 plots.

Pulled m8 via 16 chunked `read_gcs_object` calls (5 MB each) because
`download_gcs_file` was permission-denied on `/workspace`; concatenated tool-
result JSONs with `jq -j '.content'`.

**Why.**
- Thresholds anchored to score distribution sample (see
  THRESHOLD.md): bits histogram bimodal, valley near 50; fident is
  uninformative across the whole range (median ~0.18, mostly structural
  homologs at low seq id) → kept fident out of the gate.
- Host fill *layers* rather than single source: ICTV is canonical but
  family/genus only; rules + name regex + curated human override fill the
  long tail. Per HOST_FILL.md, expected unfilled ~10-15% which is what we got.
- Per-virus 14-feature representation (n_hits, mean_bits, fident stats,
  plddt, orf-source fractions) for UMAP: avoids the sparse counts-matrix
  problem that would otherwise collapse a single-sample experiment.

**Evidence.**
- Ingest: 752,332 raw → 399,013 discovery (53.0%) → 101,251 high-conf
  (13.5%). 15,511 unique queries, 24,925 unique targets, 6,586 taxids,
  185 families. Unmapped uniprot→taxid: 117 (0.03%).
- Host distribution (discovery hits): bacteria 197,107 (49.4%) ; protist
  105,089 (26.3%) ; unknown 52,371 (13.1%) ; vertebrate_nonhuman 13,547 ;
  arthropod 13,745 ; plant 9,156 ; archaea 4,348 ; fungus 2,083 ; human 1,450.
- ORF source: assembly 333,828 ; vs2 55,984 ; genomad 9,069 ; rescued 132.
  (Initial parser missed `genomad|` prefix → 9k unknowns; fixed.)
- Outputs in `out/HVP-0006.1_34P/`: anndata.h5ad (1.0 MB),
  hits_filtered.parquet, summary.json, plots/ (10 files: 8 PNG + 2 HTML).

**Bugs hit + fixes.**
- AnnData stores `var` string cols as Categorical; `df["family"].fillna(...)`
  raises `TypeError: Cannot setitem on a Categorical with a new category`.
  Fix in plots.py + plots_tier2.py: `astype("object").where(notna, ...)`
  before any string replacement.
- `gcloud` / `gsutil` not in container; `download_gcs_file` permission-denied.
  Workaround: chunked `read_gcs_object` → `jq -j .content` concat.

**Outcome.** Single-sample pipeline runs end to end in ~1 min from m8.
All 10 plots written. Headline finding for HVP-0006.1_34P (a respiratory
sample): 49% of structural hits are to bacteriophage, 26% to protist
viruses, only 0.4% to viruses with curated `human` host — a real human
respiratory sample is dominated by phage of resident bacteria and
environmentally-coincident viruses, not human viruses directly. Worth
checking on more samples before drawing any cohort-level claim.

**Follow-up.**
- Cohort scaffold: pipeline is single-sample today. Need driver that
  ingests N samples → samples×viruses matrix → real UMAP across samples.
- Foldseek `makepaddedseqdb` padding still pending — separate work stream.

## 2026-06-23T22:18:00Z — viral-viz: host-fill v2 (UniProt + class/order rules + NCLDV regex)

**Context.** v1 host fill (ICTV + family rules + name regex) left 52,371 hits
(13.1% of 399k) as `unknown` for HVP-0006.1_34P. Top blockers: NCLDV giants
(Pithovirus/Pandoravirus/Pacmanvirus/Faustovirus/Kaumoebavirus) lacking family,
"Prokaryotic dsDNA virus sp." (10k hits alone) lacking lineage, and taxa with
class/order set but family null (Pisoniviricetes/Bunyavirales/Picornavirales).

**Decision / action.** Layered v2:
- `uniprot_host.py` → batched UniProt REST (100 per query)
  against the sample's 24,159 valid UniProtKB accessions (UPI* filtered).
  Cache `refs/cache/uniprot_host.parquet`.
- `host_taxid_lookup.py` → maps cited host taxids → host_group
  via NCBI eutils efetch on `db=taxonomy` (XML). Curated `KNOWN_HOSTS` covers
  common model organisms first; rest get lineage-bucketed.
- `build_host_table.py` rewritten with 4 layers (L1 ICTV,
  L2 UniProt majority-vote per viral taxid, L3 expanded lineage rules
  class+order+family, L4 name regex). Class/order rules derived from ICTV
  MSL41 majority-vote (frac ≥ 0.7).

**Why.**
- UniProt has `virus_hosts` cross-reference (host species + taxid) for any
  characterized viral entry — but most BFVD entries (uncharacterized
  metagenomic predictions) won't have it. We found out only 480/24,159
  (2.0%) of sample uniprots actually carry host info — small fish.
- The bigger lever was class/order rules + NCLDV genus regex: most of the
  remaining unknowns *do* have class or order assigned, just not family.
- Per-viral-taxid majority vote handles uniprots that cite multiple hosts.

**Bugs hit + fixes.**
- UniProt REST 400'd on UPI* (UniParc) accessions; added UNIPROTKB_RE
  filter. Also added bisect-on-400 fallback in case other poisoned IDs
  slip through.
- `lineage_to_group()` initially keyed on `superkingdom`, returned all
  Eukaryotes as unknown. NCBI taxonomy now emits `domain` instead;
  fallback to either fixed it (`refs/cache/host_taxid_to_group.parquet`
  jumped from 26 mapped to 236 of 238).

**Evidence.** Hits-by-host_group before → after on HVP-0006.1_34P:

| host_group           | before  | after   | Δ        |
|----------------------|--------:|--------:|---------:|
| bacteria             | 197,107 | 207,543 |  +10,436 |
| protist              | 105,089 | 129,984 |  +24,895 |
| unknown              |  52,371 |  15,840 |  -36,531 |
| vertebrate_nonhuman  |  13,547 |  14,484 |     +937 |
| arthropod            |  13,745 |  13,745 |        0 |
| plant                |   9,156 |   9,419 |     +263 |
| archaea              |   4,348 |   4,348 |        0 |
| fungus               |   2,083 |   2,083 |        0 |
| human                |   1,450 |   1,450 |        0 |
| nan                  |     117 |     117 |        0 |

Unknown share dropped 13.1% → 4.0% of all discovery hits. Most of the
protist gain (+24,895) is NCLDV genera (Pithovirus/Pacmanvirus/Faustovirus/
Kaumoebavirus/etc.); bacteria gain (+10,436) is mostly "Prokaryotic dsDNA
virus sp." (taxid 2591644). Top 30 abundance bar + host UMAP both shift —
no high-confidence human hit moved (still 1,450 = 0.4%), so the headline
conclusion (sample dominated by phage/protist viruses, not human viruses)
holds.

**Outcome.**
- `refs/cache/uniprot_host.parquet`: 24,159 rows (480 with host string)
- `refs/cache/host_taxid_to_group.parquet`: 238 rows (2 unknown residual)
- `refs/taxid_host.parquet` rebuilt; host_source distribution now includes
  `uniprot_host` (5), `rule_order` (381), `rule_class` (3,642) layers.
- All 10 plots regenerated; v1 versions preserved in
  `out/HVP-0006.1_34P/_before_v2_host/` for visual diff.

**Follow-up.** Open: 15,840 still-unknown hits are dominated by
"uncultured virus" / "uncultured marine virus" entries with no lineage at
all — would need read-level metadata or a CRISPR-spacer host prediction
(iPHoP/PHIST) pass to push further. Out of scope unless residual matters
for cohort analysis.

## 2026-06-23T22:24:54Z  — UMAP fix: scalar-feature → virus×ORF sparse matrix

**Context.** User feedback on v2 plots: "Umaps look bad. Can you fix it?"
The previous `_build_virus_feature_matrix` produced 14 scalar hit-quality
features per virus (mean/median/std of bits, fident, plddt, etc.). These
carry zero taxonomic or host signal — they describe how well foldseek
scored hits, not which ORFs hit which virus. UMAP can't recover structure
that isn't in the input features.

**Decision / action.** Replaced `_build_virus_feature_matrix` →
`_build_virus_orf_matrix` in `plots.py:145`. Builds
sparse virus×ORF count matrix (single-cell analog: viruses=cells,
ORFs=genes, hits=expression). Pipeline: median library-size normalize
→ log1p → TruncatedSVD(30) → L2 norm → UMAP(cosine, n_neighbors=30,
min_dist=0.1). Plot styling: legend outside chart, unknown-host points
drawn first as faint background then colored hosts overlay.

**Why.** Two viruses are biologically similar when their shared-ORF
profiles overlap (homology signal), not when their foldseek score
distributions are similar. Sparse matrices are the right shape for
scanpy-style workflows: density 0.475% on this sample, well within
SVD+UMAP comfort zone. n_neighbors=30 / min_dist=0.1 are the scanpy
defaults — known good for biological count data.

**Evidence.**
```
[plot] building virus × ORF sparse matrix
[plot]   shape=(4892, 12479)  nnz=289,914  density=0.475%
[plot] TruncatedSVD → 30 components
[plot]   explained variance: 30.3%
[plot] UMAP (n_neighbors=30, metric=cosine)
[plot]   wrote out/HVP-0006.1_34P/plots/umap_by_family.png
[plot]   wrote out/HVP-0006.1_34P/plots/umap_by_host.png
```
Filters dropped viruses with <3 hits and ORFs hitting <2 viruses (these
add noise — singletons contribute one row with no neighbors). 4,892 of
6,587 viruses survive (74%).

**Outcome.** New PNGs at `out/HVP-0006.1_34P/plots/umap_by_{family,host}.png`.
Awaiting user review of cluster quality.

## 2026-06-23T22:30:00Z  — UMAP: add leiden clustering (scanpy canonical pipeline)

**Context.** User: "No leiden / louvain clustering?" After UMAP fix landed,
the scatter showed visible structure (handful of tight blobs + diffuse
manifold) but no unsupervised cluster labels — only family / host overlays.

**Decision / action.** Migrated `plot_umap` from raw `umap-learn` to the
canonical scanpy pipeline: build AnnData(obsm["X_pca"]=Xz) → sc.pp.neighbors
(n_neighbors=30, cosine, use_rep="X_pca") → sc.tl.leiden(resolution=1.0,
flavor="igraph", n_iterations=2) → sc.tl.umap(min_dist=0.1). Both leiden
and UMAP now reuse one kNN graph instead of building two separately.
Added third output `umap_by_cluster.png` colored by leiden cluster id;
top-20 clusters distinct colors, rest merged to "other" grey.

**Why.** Two reasons. (1) Sharing the neighbor graph between clustering and
embedding is the scanpy convention — leiden labels and UMAP positions are
then guaranteed coherent, no "cluster split across two UMAP islands"
artifact. (2) scanpy's leiden is faster and more robust than calling
leidenalg directly. resolution=1.0 is scanpy default — good middle ground
between over-splitting (many tiny clusters) and under-merging (everything
in one blob). On this sample it landed at 45 clusters across 4,892 viruses
(median ~109 viruses/cluster).

**Evidence.**
```
[plot] scanpy neighbors (n_neighbors=30, metric=cosine) + leiden + UMAP
[plot]   leiden: 45 clusters
[plot]   wrote out/HVP-0006.1_34P/plots/umap_by_cluster.png
```

**Outcome.** Three UMAP PNGs: by_family / by_host / by_cluster. Awaiting
user review of cluster–family / cluster–host correspondence (whether
unsupervised structure aligns with the curated labels, or surfaces
novel groupings worth investigating).

## 2026-06-23T22:42:00Z  — Cluster labeler + biological interpretation

**Context.** User: "Can you find some assignment to the clusters that
makes logical sense?" Need to translate 45 leiden cluster ids into
biologically meaningful names.

**Decision / action.** Wrote `label_clusters.py`.
Reads `plots/clusters.parquet` (virus × cluster + metadata), summarizes
each cluster: dominant family + purity, dominant host + purity, top-3
exemplar scientific names *within the dominant family* (fixes a
giant-virus bias — initial version picked exemplars by hit_count which
surfaced Pithovirus/Pandoravirus on every cluster they touched), and a
GIANT flag when top-10% of cluster members carry >60% of hits. Also
modified `plots.py:plot_umap` to persist `clusters.parquet` alongside
the PNGs.

**Why.** Raw "c0..c44" ids are useless for a writeup. ICTV families
are the right reference vocabulary, but most clusters are mixed at the
family level — phages especially. Purity buckets (≥50% = name it,
≥25% = "X-leaning", <25% = "polyphyletic") let the label honestly
reflect cluster coherence instead of pretending purity that isn't there.

**Evidence.** 45 clusters across 4,892 viruses. Clean stratification:
- 7 single-family clusters (purity ≥60%): Retroviridae, Parvoviridae,
  Casjensviridae, Microviridae, Straboviridae, Herelleviridae,
  Autographiviridae.
- ~17 phage super-clades (host=bacteria ≥80%, family <50%) — partition
  Caudoviricetes by shared packaging/capsid module.
- 5 eukaryotic-RNA-virus clades (Coronaviridae, Astroviridae,
  Closteroviridae, Endornaviridae, Baculoviridae).
- NCLDV + virophage clusters (c0 Poxviridae+Mimiviridae;
  c33 Lavidaviridae; c1 Phycodnaviridae+Kyanoviridae).
- Cross-domain ORF bridges (c8 +ssRNA picorna-superfamily polyprotein;
  c5/c11/c12 herpes+phage via HK97 capsid fold).

**Outcome.** `plots/cluster_labels.tsv` written, summary surfaced to
user. Key insight to flag in writeup: clusters track
**replication/packaging machinery homology**, not strict ICTV family.
The herpes-in-phage clustering is real biology (HK97 capsid fold is
conserved between tailed dsDNA phage and herpesviruses) — not a UMAP
artifact.

## 2026-06-24T00:00:00Z  — Leiden resolution scan: rank-coherent ICTV order

**Context.** User: "Optimize the number of clusters based on some
standardized measure? I don't like mixing hierarchy values." The
res=1.0 default produced 45 clusters mixing pure ICTV families
(c22 Retroviridae 96%) with phage super-clades (Caudoviricetes-spanning,
family <25%) — clusters at inconsistent hierarchy depth.

**Decision / action.** Wrote `scan_resolution.py`.
Refactored `plot_umap` to expose `compute_virus_embedding(hits)` returning
adata with kNN graph + UMAP but no leiden — scan reuses one embedding
across 11 resolutions (0.1 → 4.0). For each, computed: Newman modularity
on the connectivities graph, silhouette on the SVD embedding (cosine,
sampled 2000), and ARI + NMI vs ICTV class/order/family/genus
(restricted to viruses where the rank is annotated). Picked
rank-coherent (resolution, rank) by max ARI with NMI tiebreaker.
`label_clusters.py` gained `--rank` flag so labels match the recommended
hierarchy level; `plots.py` now persists class/order/family/genus
columns alongside cluster ids.

**Why.** Modularity rewards more clusters and is flat (0.87 → 0.93)
across the whole range — useless as a discriminator. Silhouette peaks
at res=0.7 (0.566) but tracks geometric separation, not biological
coherence. ARI is the right metric: it's normalized for cluster count
AND it directly measures alignment to a known partition (the ICTV
hierarchy). Picking the **single rank** where ARI peaks forces all
clusters to live at one hierarchy level — exactly what the user
asked for.

**Evidence.**
```
res=0.10  k= 15  mod=0.872  sil=0.338  ARI(class/order/family/genus)=0.065/0.303/0.099/0.008
res=0.30  k= 32  mod=0.923  sil=0.513  ARI=0.031/0.376/0.154/0.022   ← PICK
res=0.70  k= 42  mod=0.926  sil=0.566  ARI=0.028/0.375/0.160/0.024
res=1.00  k= 47  mod=0.928  sil=0.557  ARI=0.026/0.284/0.170/0.024
res=4.00  k= 85  mod=0.909  sil=0.448  ARI=0.021/0.270/0.170/0.033
```
Order ARI peaks at res=0.3 (0.376, NMI=0.560). Class ARI is uniformly
tiny (~0.03) — too few classes to discriminate; this is a healthy
sanity check (clusters aren't just sitting at the top of the tree).
Genus ARI is uniformly tiny (~0.02) — genus is too fine for ORF-sharing
signal at this depth.

**Outcome.** Wrote `plots/leiden_scan.{tsv,png}` and regenerated
`umap_by_*` PNGs at res=0.3 with order-rank labels. 32 clusters; the
pure-order standouts: c19 Ortervirales/vertebrate (100%/87%, HIV+SIV+HTLV),
c27 Lefavirales (100%), c28 Algavirales (100%), c29 Kirjokansivirales
(100%), c1 Chitovirales (100%), c2 Piccovirales (70% — parvoviruses),
c11 Martellivirales/plant (52% — leafroll viruses), c21 Picornavirales
(68%). The previous "phage super-clades" at family rank now resolve to
named Caudoviricetes orders (Crassvirales, Petitvirales, Tubulavirales,
Kirjokansivirales, Algavirales, Lefavirales, Methanobavirales) — the
ICTV reorganization of phage taxonomy into orders matches what leiden
finds in the ORF-sharing graph.

## 2026-06-24T12:20:00Z  — Protein-level cluster annotation, opt-in flag

**Context.** Previous session shipped leiden clustering + ICTV-rank cluster
labels for `HVP-0006.1_34P`. User asked: "Verify that you're getting these
labels from the protein information." Answer: the **signal** is protein-level
(foldseek 3Di structural-similarity hits per query ORF) but the **labels**
are derived from taxid joins (host_group, ICTV order). Actual protein
function names (Pfam / GO / protein_name) were not yet wired in. User then
said: "Yes. Add them but keep a flag so we can recreate the plots and
analysis as it currently is."

**Decision / action.**
1. New module `uniprot_protein_name.py` — mirrors
   `uniprot_host.py` plumbing (batched UniProt REST, `UNIPROTKB_RE` filter,
   bisect-on-400, incremental cache). Fetches `fields="accession,protein_name"`.
   Caches to `refs/cache/uniprot_protein_name.parquet`.
2. New module `cluster_markers.py` — loads
   `plots/clusters.parquet` + `hits_filtered.parquet`, rebuilds virus×ORF
   matrix via `plots._build_virus_orf_matrix`, runs
   `sc.tl.rank_genes_groups(method="wilcoxon")` per leiden cluster, joins each
   marker query to its best foldseek hit and the corresponding protein name.
   Walks down the hit list when the top-bits hit returns
   `protein_name="deleted"` (TrEMBL obsolete entries — ~70% of top hits in
   this sample). Emits `cluster_markers.tsv` (one row per marker) and
   `cluster_markers_summary.tsv` (one row per cluster, top-5 names).
3. `label_clusters.py`: added `--with-protein-markers` flag (default OFF
   preserves prior behavior). When set, calls the marker pipeline and renders
   `umap_by_cluster_labeled_<rank>_with_proteins.png` alongside the original
   PNG without overwriting it.

**Why.** The user wanted a quality check on the cluster labels using a
*different* data axis than the one that produced them. Taxonomic labels can
be wrong/noisy where the BFVD taxid→ICTV joins are incomplete; protein names
are an independent ground truth. Keeping it behind a flag preserves the
ability to regenerate prior figures byte-for-byte and isolates the
opt-in dependency on the protein-name cache.

**Evidence.**
- UniProt fetch: 24,159 accessions in 117 s (~205/s), 24,159 names cached.
- `cluster_markers.tsv`: 320 rows (32 clusters × 10 markers).
- Cluster ↔ protein content sanity:
  - c16 Herpesvirales/bacteria → Terminase large subunit (canonical tailed
    phage packaging gene; confirms the cluster is phages, not real
    herpesviruses).
  - c19 Ortervirales/vertebrate_nonhuman → Gag-Pro-Pol polyprotein (correct
    retrovirus signature).
  - c12 Tubulavirales/bacteria → RecA-like DNA recombinase.
  - c9 Herpesvirales-leaning/bacteria → Ribonucleoside-diphosphate reductase
    (phage replication enzyme).
  - c27 Lefavirales/bacteria → Integrase (lysogeny).
  - c21 Picornavirales/mixed-host → Protease Do-like PDZ + Vgr OB-fold
    (suggests mis-grouped Type VI secretion / phage-like protease, not
    picornavirus).
  - c11 Marteilevirales/plant → DNA 3'-5' helicase + Replicase polyprotein
    1ab (ORF1ab is a coronavirus replicase — mis-labeled as Marteilevirales).

**Outcome.** Protein-level annotation closes the loop between embedding
signal and label. Cluster labels stay taxonomy-only by default; users get
protein content as an opt-in overlay. New artifacts in
`out/HVP-0006.1_34P/plots/`:
- `cluster_markers.tsv`, `cluster_markers_summary.tsv`,
  `umap_by_cluster_labeled_order_with_proteins.png`.

Original artifacts (`cluster_labels_order.tsv`,
`umap_by_cluster_labeled_order.png`) untouched.

**Follow-up.** Several clusters (c11 "Marteilevirales", c21 "Picornavirales")
show protein content that contradicts the dominant ICTV order — the BFVD
taxid→ICTV-order join must be sparse for those targets. Worth verifying
whether the protein names are the correct label (i.e. patch ICTV joins) or
whether these proteins are real cross-host homologs. Not blocking — current
pipeline now has both signals visible side-by-side.
