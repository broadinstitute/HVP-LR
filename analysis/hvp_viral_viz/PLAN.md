# Implementation plan — viral-viz

Sequenced. Each phase produces a concrete artifact + check.

## Phase 0 — Foundations (done as part of this plan write-up)

- [x] Folder skeleton + docs (README, PLAN, THRESHOLD, HOST_FILL)
- [x] Confirm BFVD has taxonomy lineage with realm→species (351,125 rows)
- [x] Confirm host species NOT in BFVD; external fill needed
- [x] Empirical sample of foldseek m8 score distribution
- [x] Pick threshold + write rationale (see THRESHOLD.md)

## Phase 1 — Reference tables (refs/)

1. **`build_bfvd_refs.py`** — collapse BFVD into per-uniprot tables.
   Inputs: 3 BFVD TSVs from `/workspace/bfvd/`.
   Outputs in `refs/`:
   - `uniprot_taxid.parquet` — uniprot_id, taxid (≈ 347,514 rows)
   - `taxid_lineage.parquet` — taxid, scientific_name, realm, kingdom, phylum, class, order, family, genus, species (≈ 23,614 rows; parse `d_/k_/p_/c_/o_/f_/g_/s_` prefixes)
   - `uniprot_qc.parquet` — uniprot_id, model_id, avg_plddt, ptm, splitted, version (best per uniprot when split)
   Validate: every BFVD uniprot has taxid; ≥ 99% of taxids resolve to realm + species.

2. **`build_host_table.py`** — produce taxid→host_group.
   Strategy in `HOST_FILL.md`. Outputs:
   - `taxid_host.parquet` — taxid, host_group, host_source (rule|ictv|uniprot|ncbi|unknown), host_organism
   - `host_fill_audit.tsv` — per-taxid which source filled it

## Phase 2 — Hit ingestor (src/hvp_viral_viz/ingest.py)

3. **`ingest_run`** — CLI consuming one m8 + sample metadata:
   - Parse m8 (7 cols: query, target, evalue, bits, fident, alnlen, mismatch)
   - Apply discovery threshold (THRESHOLD.md → bits ≥ 50 AND evalue ≤ 1e-5)
   - Join `target` → uniprot_id → taxid → lineage → host
   - Tag each query with ORF source from query prefix (`vs2|`, `assembly|`, `rescued|`)
   - Emit `out/per_sample/<sample>.hits.parquet`

4. **`build_anndata`** — combine N per-sample parquets:
   - **Framing A (cohort, primary):** virus × sample matrix
     - obs = unique target uniprots (or collapsed to taxid for sparser viz)
     - var = sample_id (current = 1 sample)
     - .X = sum bits (or count hits, or binary) per (target, sample)
     - .obs columns: taxid, scientific_name, family, genus, host_group, avg_plddt, ptm
     - .var columns: sample_id, entity_type, cohort, donor, n_queries, n_orfs_vs2, n_orfs_assembly, n_orfs_rescued
   - **Framing B (single-sample fallback):** virus × query matrix
     - obs = unique targets, var = unique queries, .X = bits
     - Only computed for samples flagged as singleton viz target
   - Save as `out/anndata/cohort.h5ad` + `out/anndata/<sample>.singleton.h5ad`

## Phase 3 — Tier 1 plots (user-requested)

5. **Count matrix** — `plot_count_heatmap.py`
   - Clustered heatmap (viruses × samples), rows=top-N by total bits
   - CSV companion: `out/tables/virus_x_sample_counts.csv`

6. **UMAP × taxonomy** — `plot_umap_taxonomy.py`
   - Cohort mode (Framing A): `sc.pp.normalize` → `sc.pp.pca` → `sc.pp.neighbors` → `sc.tl.umap`; color by family / genus
   - Single mode (Framing B): same pipeline on virus×query matrix
   - One PNG per taxonomy rank: realm, phylum, class, family, genus

7. **UMAP × host** — `plot_umap_host.py`
   - Reuse coords from #6; recolor by host_group from `taxid_host`

8. **Abundance bar** — `plot_abundance_bar.py`
   - Top-N (default 30) viruses by hit count, stacked by family
   - Also: top-N by summed bits

## Phase 4 — Tier 2 diagnostic + interpretive plots

9. **Krona / sunburst** — `plot_taxonomy_sunburst.py` (plotly)
10. **ORF-source Sankey** — `plot_orf_source_sankey.py` (vs2/assembly/rescued → cluster → family)
11. **bits×fident scatter** — `plot_quality_scatter.py` colored by family
12. **plDDT distribution per family** — `plot_plddt_per_family.py`
13. **Hit-coverage histogram** — `plot_coverage_histogram.py` (alnlen / query_len)

## Phase 5 — Tier 3 biology-first (deferred unless asked)

14. Per-genome reconstruction grid — only when ≥ 1 target has many query hits
15. Phylogenetic placement — out of scope for v1
16. Bipartite network — out of scope for v1

## Phase 6 — CLI + reproducibility

17. **`hvp-viral-viz` console entry** — `pip install -e .` exposes:
    ```
    hvp-viral-viz refs build
    hvp-viral-viz ingest <m8> --sample <id> --entity-type <type>
    hvp-viral-viz anndata build --out cohort.h5ad
    hvp-viral-viz plot all --anndata cohort.h5ad --out out/
    ```
18. **`Makefile`** at viral-viz root with `make refs`, `make ingest`, `make plots`, `make all`
19. **`environment.yml`** — scanpy, anndata, plotly, pyarrow, requests
20. **`tests/`** — unit tests on small synthetic m8 + small BFVD subset

## Phase 7 — Future cohort

21. Re-run `ingest_run` per new sample; `build_anndata` rebuilds cohort h5ad
22. Sample-level PCoA + diversity plots gated on n_samples ≥ 3
23. Per-virus prevalence vs mean-bits plot gated on n_samples ≥ 5

## Open questions to close before coding Phase 2

- best-hit per (query, target) or sum all (query, target) bits? → sum bits in `.X`, also store binary 1/0 layer
- count "virus" at uniprot, taxid, or species level? → keep uniprot as primary obs key, derive species-level rollup as secondary AnnData
- treat splitted models (`A0A2P1GMZ4_1` vs `A0A2P1GMZ4`) as same virus? → yes; collapse on uniprot stem
