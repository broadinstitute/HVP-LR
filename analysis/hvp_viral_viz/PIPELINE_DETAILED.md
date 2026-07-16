# HVP-LR cohort viral landscape — pipeline (detailed)

How PacBio HiFi sequencing reads become the UMAP scatter plots of the virome.

The cohort is 48 specimens (24 saliva + 24 stool) from the Human Virome Project long-read arm. Each specimen is processed independently through a series of WDL pipelines on Terra; per-specimen results are pooled and visualized by the local `hvp_viral_viz` package. The end-to-end workflow chain is:

```
HvpReadProcessing → HvpAssembly → HvpReadRescue → HvpMagsPipeline → HvpViralPipeline → HvpViralProteinAnnotation
```

`HvpMagsPipeline` produces the bacterial MAG catalog (parallel sibling to the viral track — its outputs do not feed the viral pipeline directly, but it is part of the per-specimen processing chain). The two workflows that produced the data feeding the UMAP plots are **`HvpViralPipeline`** and **`HvpViralProteinAnnotation`**; everything else is per-specimen context.

## A. Upstream — from sequencer to merged viral input FASTA

1. **PacBio HiFi sequencing.** Each specimen is sequenced on a PacBio long-read instrument with circular-consensus (CCS) chemistry. Output is a per-sample unaligned HiFi BAM of Q20+ reads in the 10–25 kb range.

2. **`HvpReadProcessing`.** HiFi reads are cleaned of PhiX/spike-ins and human host reads, then classified with **Kraken2** for per-read taxonomic status.

3. **`HvpAssembly`.** The cleaned reads are assembled with **Myloasm**, a HiFi-tuned long-read metagenome assembler. Reads are also mapped back to the assembly with **minimap2** (`map-hifi`) producing a coordinate-sorted BAM.

4. **`HvpReadRescue`.** The Kraken2 per-read classifications are cross-referenced against the minimap2 read-to-contig BAM. Reads that are *not* accounted for by the assembly — unmapped or with soft+hard clip fraction ≥ 0.5 — are partitioned into a **kraken-viral rescue set** and a **kraken-unclassified rescue set** (default minimum rescue-read length: 1 kb). The assembly contigs plus both rescue sets are concatenated into a single **merged FASTA** (`merged_fa_gz`). This is the primary nucleotide input to `HvpViralPipeline`.

5. **`HvpMagsPipeline` (parallel bacterial track).** On the same primary assembly + sorted BAM from `HvpAssembly`, three binners run in parallel — **MetaBAT2**, **MaxBin2**, **SemiBin2** (env-specific model: `human_gut` for stool, `human_oral` for saliva) — sharing a single per-contig depth profile from `JgiDepth`. **DAS_Tool** dereplicates and refines the union of bins (score threshold 0.6). Each refined bin is then evaluated by **CheckM2** (completeness, contamination) and assigned a GTDB r226 taxonomy via **skani**. Per-sample MAG scalars (HQ/MQ/LQ bin counts, distinct GTDB species/genera, % assembly bases binned) are emitted for the Terra data table. This is the **bacterial host catalog**; it is independent of the viral pipeline (the viral track uses `merged_fa_gz` from `HvpReadRescue`, not the MAG bins) but provides the per-specimen microbial context.

## B. `HvpViralPipeline` — per-specimen viral contig identification

6. **t_01: `Genomad`.** The merged FASTA is run through geNomad's neural-network viral-classifier. Outputs: `virus_summary.tsv` (per-sequence scores, topology, taxonomy), `virus.fna` (viral nucleotide FASTA), **`virus_proteins.faa` (geNomad's own predicted protein FASTA — passed straight to the protein-annotation pipeline without an additional ORF call)**, and a per-gene TSV.

7. **t_02: `VirSorter2`.** Runs in parallel with geNomad on the same merged FASTA. Outputs: `final-viral-score.tsv` (per-sequence viral scores and group assignments) and `viral_combined.fa` (viral nucleotide FASTA). Default `min_length` = 1000 bp.

8. **t_03, t_04: `CheckV` (× 2).** CheckV runs once per upstream tool's nucleotide output (geNomad and VS2 separately), producing `quality_summary.tsv` plus separate `viruses.fna` (complete / high-quality genomes) and `proviruses.fna` (extracted proviral regions).

9. **t_05: `ViralContigSummary`.** Joins all four upstream outputs (geNomad summary + VS2 score + both CheckV summaries) into a single per-contig TSV indexed by contig ID.

10. **t_06: `ViralOverallSummary`.** Aggregates the per-contig table to sample-level scalar metrics (number of viral contigs, number of HQ/complete/MQ/LQ genomes, N50 of viral contigs, total viral bases, mean completeness/contamination, DTR/ITR/provirus counts).

## C. `HvpViralProteinAnnotation` — per-specimen structural protein search

Inputs: VS2 `viral_combined.fa` (always), geNomad `virus_proteins.faa` (always; pre-computed by geNomad), optionally an assembly-contigs FASTA, optionally rescued-reads FASTA. Plus BFVD foldseek DB (+ padded variant for GPU), ProstT5 weights, BFVD reference-metadata TSV.

11. **t_01–t_03: `PyrodigalGv` × up to 3.** Pyrodigal-gv (viral-genome-aware Prodigal fork) calls ORFs on the VS2 viral combined FASTA (always), and optionally on the assembly contigs and the rescued reads if those inputs are supplied. The geNomad pre-computed protein FASTA is concatenated directly — **no additional ORF call**.

12. **t_04: `ConcatProteinFastas`.** All AA FASTAs are concatenated into one, with **per-source header prefixes** — `genomad|`, `vs2|`, `assembly|`, `rescued|` — so every downstream hit remains attributable to its ORF source.

13. **t_05: `MmseqsEasyLinclust` (NR collapse).** Concatenated proteins are clustered with `mmseqs easy-linclust` at **90% AA identity, 80% coverage** (defaults). One representative sequence per cluster is kept; this is the **NR protein set** that proceeds to the structural search. Cluster-membership TSV is preserved so each NR representative can be expanded back to all source ORFs.

14. **t_06: `FoldseekCreateDbFromFasta`.** The NR protein FASTA is converted into a foldseek structure DB via **ProstT5** (a protein-language model that translates each amino-acid sequence into a 1-D structural alphabet — 3Di tokens). No physical or AlphaFold structure is required. GPU pathway (default) runs on `foldseek-gpu` images on `nvidia-tesla-t4`.

15. **t_07: `FoldseekSearch` vs BFVD.** Structural search of the NR protein DB against **BFVD** (Big Fantastic Virus Database — ~350k AlphaFold-predicted viral protein structures, indexed by UniProt accession). GPU search requires the **padded BFVD DB** (produced once via `foldseek 10 makepaddedseqdb`); CPU search uses the unpadded DB. Default `evalue_cutoff` = 0.001.

16. **t_08: `FoldseekConvertAlis`.** Alignment DB is formatted as a TSV with columns `query, target, evalue, bits, fident, alnlen, mismatch`. (Note: queries here are ProstT5-built sequence-only foldseek DBs — they have no CA datafile, so structural-coordinate columns are unavailable.)

17. **t_09: `AnnotationTransfer`.** Best hit per query (lowest e-value, ties broken by higher bit score). Left-joined on `target` against the BFVD reference-metadata TSV — columns: `uniprot_id, model_id, avg_plddt, ptm, splitted, version`. Output: `best_hits_tsv` (one row per query) and `annotated_hits_tsv` (best hit + reference metadata).

## D. Post-pipeline local analysis (`hvp_viral_viz`)

18. **Ingest hits (`ingest.py`).** Per-sample `foldseek_hits.tsv` is read, then filtered to a **discovery tier**: `evalue ≤ 1e-5 AND bits ≥ 50 AND alnlen ≥ 50`. A **high-confidence tier** (`bits ≥ 300 AND alnlen ≥ 80 AND evalue ≤ 1e-10`) is also computed and used for taxonomic anchors. Output: `hits_filtered.parquet`.

19. **BFVD → UniProt → taxid → ICTV lineage.** Each foldseek `target` is mapped to its UniProt accession (the BFVD metadata `uniprot_id` column), then to a host taxid, then joined against the ICTV Virus Metadata Resource (VMR) to attach class, order, family, genus, species. The ICTV **order** is the primary grouping rank used downstream.

20. **Host-group fill — five-tier cascade (`build_host_table.py`).** Each taxid's host group (`bacteria`, `archaea`, `vertebrate_human`, `vertebrate_nonhuman`, `invertebrate`, `plant`, `protist`, `fungi`, `unknown`) is filled by walking these sources in order:
    1. ICTV VMR `host_range` field.
    2. UniProt `organism_host` metadata.
    3. NCBI Virus host lookup.
    4. Rule-based heuristics on the lineage (e.g. *Caudoviricetes* → bacteria).
    5. Virus-name pattern matching (e.g. `*phage*` → bacteria, `*-pox virus` → vertebrate).
    First source with a non-null answer wins.

21. **Virus × ORF count matrix (`plots._build_virus_orf_matrix`).** The filtered hits are pivoted into a sparse integer matrix `X` of shape *(viruses × ORFs)*, where `X[v, o]` is the count of ORFs from viral contig `v` whose best Foldseek hit is to BFVD reference protein `o`. Direct single-cell-genomics analogy: viruses are "cells", reference proteins are "genes".

22. **Normalize.** Total-count normalization to fixed library size, then `log1p` — standard scanpy preprocessing.

23. **SVD → cosine kNN → UMAP.** Truncated randomized SVD to 50 components; cosine k-nearest-neighbour graph (k = 15); UMAP embedding (default 2-D parameters).

24. **Leiden community detection (`scan_resolution.py`).** Cluster the kNN graph with the Leiden algorithm. **Resolution is selected against ICTV taxonomy**: scan γ from 0.1 to 4.0, score each clustering by adjusted Rand index (ARI) against every ICTV rank (class, order, family, genus), pick the (γ, rank) pair maximizing ARI; NMI is the tiebreaker. The chosen resolution is the one where data-driven clustering and ICTV taxonomy align best.

25. **Cluster labeling (`label_clusters.py`).** Each cluster gets a human-readable label `<dominant_order> / <dominant_host_group>` with a purity tag:
    - `≥ 50%` dominant → bare label (`Caudoviricetes / bacteria`).
    - `25–50%` → `-leaning` suffix (`Chitovirales-leaning / vertebrate_nonhuman`).
    - `< 25%` → `polyphyletic / <host>`.

26. **Marker ORFs (`cluster_markers.py`).** For each cluster, `scanpy.tl.rank_genes_groups` (Wilcoxon rank-sum) finds the top-K ORFs most enriched in that cluster vs all others. Each marker ORF's best Foldseek hit is joined to UniProt's `protein_name` — fallback walks down the hit list when the top-bits hit returns `deleted` / `Uncharacterized protein`. Outputs the named viral proteins (capsid, polymerase, …) characterizing each cluster.

## E. Cohort-level views

27. **Pool (`pool_cohort.py`).** All per-specimen `hits_filtered.parquet` files are concatenated into a single `pooled_hits.parquet` (~93M rows for the 48-specimen cohort). The same SVD → kNN → UMAP → Leiden pipeline is rerun on the pooled matrix to produce a cohort-level embedding.

28. **Sub-cohort splits.** The pool is rebuilt for the **saliva-only (24)** and **stool-only (24)** sub-cohorts. Each has its own UMAP, Leiden clusters, and labels — allowing direct comparison of which viral families/orders dominate each biospecimen type.

## F. Supplementary plots

29. `umap_by_cluster_labeled_order.png` — cluster scatter labeled with order/host strings.
30. `umap_by_family.png`, `umap_by_host.png`, `umap_by_prevalence.png` — by-rank colourings.
31. `heatmap_orf_source_top30.png` — per-cluster ORF-source breakdown (vs2 / genomad / assembly / rescued), exploiting the source-label prefix added in t_04 of `HvpViralProteinAnnotation`.
32. `abundance_bar_top30.png`, `abundance_by_family_top30.png` — top-30 prevalence bars.
33. `sample_diag_*` — per-specimen PCA, hierarchical clustering of specimens, by-specimen / by-specimen-type UMAP.

---

## Tools and references

- **Sequencing.** PacBio HiFi (Revio / Sequel IIe).
- **Read processing / assembly.** Kraken2; Myloasm (HiFi long-read metagenome assembler); minimap2 `map-hifi`.
- **Read rescue.** Custom WDL task — Kraken2 + BAM cross-reference (`HvpReadRescue`).
- **Binning / MAG QC (parallel bacterial track).** [MetaBAT2](https://bitbucket.org/berkeleylab/metabat), [MaxBin2](https://sourceforge.net/projects/maxbin2/), [SemiBin2](https://github.com/BigDataBiology/SemiBin), [DAS_Tool](https://github.com/cmks/DAS_Tool), [CheckM2](https://github.com/chklovski/CheckM2), [skani](https://github.com/bluenote-1577/skani) + GTDB r226.
- **Viral calling.** [geNomad](https://github.com/apcamargo/genomad), [VirSorter2](https://github.com/jiarong/VirSorter2), [CheckV](https://bitbucket.org/berkeleylab/checkv).
- **ORF calling / NR collapse.** [Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv); [MMseqs2](https://github.com/soedinglab/MMseqs2) `easy-linclust` (90% AA-ID, 80% coverage).
- **Structural search.** [ProstT5](https://github.com/mheinzinger/ProstT5) (AA → 3Di), [Foldseek](https://github.com/steineggerlab/foldseek), [BFVD](https://www.bfvd.foldseek.com/).
- **Taxonomy / host.** ICTV VMR, UniProt, NCBI Virus.
- **Embedding / clustering.** scanpy (PCA / kNN / UMAP); python-leidenalg; sklearn metrics.
