# HVP-LR cohort viral landscape — pipeline overview

What the UMAP plots show, and how we got from the sequencer to them.

The cohort is 48 specimens (24 saliva + 24 stool) from the Human Virome Project long-read arm. Each scatter point in the UMAP is a single viral contig; the spatial position is a structural-similarity embedding of all the contig's proteins. Clusters group viruses that share a structural-protein profile. Each specimen is processed through a chained set of WDL workflows on Terra — read processing → assembly → read rescue → MAG (bacterial) catalog → viral identification → structural protein annotation — and the local `hvp_viral_viz` Python package then ingests the foldseek hits and produces the embedding and labels.

## 1. Sequencing

1. **PacBio HiFi long-read sequencing.** Each specimen is sequenced with PacBio circular-consensus (CCS) chemistry, producing accurate (Q20+) long reads in the 10–25 kb range. Long reads recover full viral genomes that short-read assemblers fragment.

## 2. From reads to a per-specimen viral input FASTA

2. **Read processing.** Reads are cleaned of spike-ins and human host material; every read is also classified by Kraken2 for an independent taxonomic call.

3. **Long-read metagenome assembly.** A HiFi-tuned long-read assembler (Myloasm) produces contigs without binning. Reads are mapped back to the assembly to identify which reads each contig "explains".

4. **Read rescue.** Reads that are *not* explained by the assembly — unmapped, or with heavy soft/hard clipping — are recovered. Those that Kraken2 has flagged as viral, or that Kraken2 cannot classify at all, are kept as a per-specimen rescue set and added to the assembly contigs. The result is a single FASTA per specimen: assembled viral genomes + rescued long reads of likely-viral origin. This merged FASTA is the input to the viral-identification step.

## 3. Parallel bacterial track (MAG catalog)

5. **Bacterial metagenome-assembled genomes.** On the same primary assembly, three binners (MetaBAT2, MaxBin2, SemiBin2) run in parallel; DAS_Tool dereplicates and refines the union. Each refined bin is then assessed for completeness/contamination (CheckM2) and given a GTDB r226 taxonomic assignment (skani). This is the **bacterial host catalog** for each specimen — not directly used by the viral UMAP, but it is the per-specimen microbial context against which the viral findings are interpreted.

## 4. Viral contig identification

6. **Two-tool viral calling.** The merged FASTA is run through two state-of-the-art viral classifiers (geNomad and VirSorter2) in parallel. Each tool independently flags contigs as viral and emits a viral-only nucleotide FASTA.

7. **Quality assessment.** A third tool (CheckV) estimates the completeness and contamination of each predicted viral contig from both upstream tools, separating likely-complete genomes from low-quality fragments. The per-tool calls and CheckV verdicts are joined into a per-contig viral summary.

## 5. Proteins → structural search → annotation

8. **Predict proteins.** Every viral contig has its ORFs (protein-coding genes) called with a viral-aware Prodigal variant (Pyrodigal-gv). Optional inputs — assembly contigs as a whole, or only the rescued reads — can be ORF-called separately for novel-virus discovery. geNomad supplies its own predicted protein FASTA, which is concatenated directly.

9. **Tag every protein by source.** Each protein header is prefixed with where the ORF came from (geNomad, VirSorter2, assembly, rescued reads), so every downstream hit remains attributable.

10. **Non-redundant collapse.** All proteins from all sources are clustered with MMseqs2 at 90% amino-acid identity (80% coverage). One representative is kept per cluster — this avoids re-searching near-identical proteins repeatedly.

11. **Translate proteins to "structure tokens".** Each representative protein is passed through ProstT5, a protein-language model that converts the amino-acid sequence into a 1-D structural alphabet (3Di) — a representation of what the folded 3-D structure looks like, derived directly from sequence. No experimental or AlphaFold structure is required.

12. **Search vs a reference structural protein database.** Each representative protein is searched against BFVD — a public database of ~350,000 AlphaFold-predicted viral protein structures — using Foldseek. The output is, for every protein, the closest structural neighbors in the reference, with confidence scores (e-value, bit score, identity, alignment length).

13. **Best-hit annotation.** The lowest-e-value hit per query is left-joined to BFVD reference metadata (UniProt ID, average pLDDT, pTM). This is the foundation for everything downstream.

## 6. Taxonomy and host annotation

14. **Lineage and host attached.** Every reference protein's source virus has a known taxid, so every hit carries the full ICTV lineage (class, order, family, genus, species). The host (bacteria, vertebrate, plant, etc.) is filled in from ICTV, UniProt, NCBI Virus, and rule-based fallbacks — whichever has data first.

## 7. The UMAP embedding

15. **Filter to confident hits.** The local analysis applies bit-score and alignment-length thresholds to keep hits that meaningfully imply protein-level homology.

16. **Build a virus × ORF matrix.** For each viral contig, count how many of its proteins hit each reference protein. This is structurally analogous to a single-cell gene-expression matrix — viruses are the "cells", reference proteins are the "genes".

17. **Reduce, embed, cluster.** Standard single-cell-style preprocessing (normalize, log-transform, principal-components reduction) is followed by a k-nearest-neighbour graph in protein space, a 2-D UMAP projection, and Leiden community detection in the kNN graph.

18. **Choose clustering resolution against ICTV.** Rather than picking a clustering granularity arbitrarily, the pipeline scans many resolutions and picks the one whose clusters best agree with ICTV taxonomy (measured by adjusted Rand index against every ICTV rank). The chosen resolution is the one where biology and data align best.

## 8. What the plots show

19. **Each point is one viral contig.** Position reflects protein-content similarity. Clusters group viruses with shared structural-protein profiles — typically corresponding to ICTV orders or families.

20. **Cluster labels are biological.** Each cluster is labeled with its dominant ICTV order and dominant host group, with purity flags: bare label for a clean ≥ 50% cluster, `-leaning` for 25–50%, `polyphyletic` for diffuse clusters. This makes the plot interpretable at a glance — e.g. "*Caudoviricetes* / bacteria" for a phage cluster, "*Chitovirales* / vertebrate" for a poxvirus cluster.

21. **Marker proteins.** For each cluster, a statistical test identifies the proteins most enriched in that cluster vs all others. These markers are joined to UniProt protein names — the result is a short list of named viral proteins (capsid, polymerase, integrase, …) that characterize each cluster. This is the biological "what is this group of viruses doing" answer.

## 9. Cohort-level and sub-cohort views

22. **Cohort pool.** All 48 specimens are combined to produce a single cohort-level UMAP, revealing structure that no individual specimen has resolution to see.

23. **Saliva vs stool sub-cohorts.** The pool is also rebuilt separately for the 24 saliva and 24 stool specimens, allowing direct visual comparison of which viral families dominate each biospecimen type and which clusters are shared.

---

## In one sentence

We sequence each specimen long-read, recover viral genomes from both assembly and unaccounted reads, predict every viral protein, convert each protein to a 3-D-structure-aware representation, search a large database of known viral protein structures, then cluster the recovered viruses by their shared protein profile — producing UMAP plots in which each point is a viral contig and each cluster groups structurally-related viruses by ICTV taxonomy and host.
