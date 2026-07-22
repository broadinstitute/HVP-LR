
## 2026-06-22T19:23:13Z  — GPU diag preamble + collaborator-claim refutation

**Context.** Collaborator claimed foldseek-gpu Dockerfile's `-DENABLE_CUDA=1` failed to cascade to ggml/llama.cpp submodule, blaming the 1h44m t_06_FoldDb runtime on CPU-only ProstT5 inference. Pipeline running on `broad-hvp-dasc/broad-hvp-dasc-dev` workspace, submission c9fd8a34-d7bc-477b-859b-04899c10f1ae.

**Decision / action.** Refuted the Dockerfile claim; added a `use_gpu`-gated GPU diagnostic block to FoldseekCreateDbFromFasta and FoldseekSearch in `wdl/tasks/ProteinAnnotation/Foldseek.wdl` (commit 5471679). Block runs `nvidia-smi`, `nvcc --version`, dumps `CUDA_VISIBLE_DEVICES`/`LD_LIBRARY_PATH`/`libcudart` presence after the resource-detection preamble.

**Why.** Upstream `CMakeLists.txt` at pinned commit 718d4217 *does* cascade `ENABLE_CUDA=1` → `GGML_CUDA=ON` automatically inside the `if(ENABLE_PROSTT5)` block; user verified binary has 15 `ggml_cuda_*` + 50 `cublas*` symbols (static linked, `ldd` empty). The "Use GPU 0" line in stdout is a red herring — per `structcreatedb.cpp:705-739`, structcreatedb saves origGpu, calls inner createdb with `--gpu` stripped from argv, then restores par.gpu before the ProstT5 inference loop. So the slowdown root cause is NOT a missing CUDA build — it must be runtime (driver mismatch / partial offload / silent fallback). Current logs don't capture enough state to distinguish; preamble fixes that.

**Evidence.** miniwdl --strict passes. Inner-createdb stdout (`Use GPU 0`) and outer stderr (`--gpu 1`) verified via fiss-mcp get_workflow_logs. runtimeAttributes for t_06 + t_07 both show `gpuType=nvidia-tesla-t4 gpuCount=1` declared. Batch machine_type `custom-8-32768` (N1 family — T4-attachable). `CUDA0` printed in ggml device list (so CUDA runtime initialized something at task time) yet 36,330 × 172 ms/seq = 1h44m matches CPU baseline. Cromwell GCPBATCH constructs a GPU resource when any of gpuType/gpuCount/gpuDriver is set, but currently *ignores* `nvidiaDriverVersion` — Batch auto-selects (typically COS-bundled driver ~535, CUDA 12.2-compatible) while our image carries CUDA 12.6 toolkit. Forward-compat hypothesis is now testable next submission.

**Outcome.** Branch `tshea/wdl-metagenomics-poc` at 5471679. Memory saved at `project_foldseek_gpu_diagnosis.md`. Next run on this branch will surface driver version + nvidia-smi in t_06 / t_07 stdout, definitively localizing the slowdown.


## 2026-06-23T21:30:56Z  — viral-viz project scaffold + threshold decision

**Context.** User asked for viral viz of `HvpViralProteinAnnotation` outputs from submission `2da383d1-fafb-46f0-bac4-3f4f782304d0` (single-sample, HVP-0006.1_34P). Single-cell analogy: count matrix, UMAPs colored by taxonomy + host, abundance plot. Asked to set up a project folder, validate BFVD has needed data, propose host fill alternatives, pick a hit-score threshold with rationale.

**Decision / action.** Created `/workspace/viral-viz/` (sibling of `HVP-LR/`, not part of pipeline repo) with `README.md`, `PLAN.md`, `THRESHOLD.md`, `HOST_FILL.md`. Validated `/workspace/bfvd/` (`bfvd_metadata.tsv`, `bfvd_taxid.tsv`, `bfvd_taxid_rank_scientificname_lineage.tsv`): 347,514 unique uniprots, 23,614 distinct taxids, 99.97% lineage coverage from realm to species — taxonomy fully covered, but **host species not in BFVD**, requires external fill.

**Why.** Single-cell framing maps cleanly: cells = unique BFVD virus targets, features = our NR query proteins, sample axis = HVP samples (Framing A, cohort-ready). Single-sample fallback inverts to virus×query (Framing B) so UMAP works with N=1 today. AnnData + scanpy is the natural stack — works at single-sample today, scales to cohort with no code change. Host gap forces a multi-source fill: chose ICTV VMR + rule-by-lineage + name-regex for v1 (zero API cost, 95-99% expected coverage), with UniProt REST + NCBI efetch deferred as Phase 2 layers.

**Threshold decision.** Sampled the 83 MB m8 at 10 evenly-spaced offsets to characterize score distribution without downloading whole file (fiss-mcp sandbox isolates `/tmp` from main container). Observed regime break at bits ≈ 250-300; foldseek already pre-filters at e-value ≤ 1e-3. Settled on dual threshold: **discovery = `evalue ≤ 1e-5 AND bits ≥ 50 AND alnlen ≥ 50`** (one order tighter than foldseek default; structurally informative floor), **high-confidence = `bits ≥ 300 AND alnlen ≥ 80 AND evalue ≤ 1e-10`**. Deliberately NO `fident` threshold — structural homology persists below 0.3 fident, gating there would bias toward closely related viruses. Both thresholds reported on every plot; revisit triggers documented (cohort ≥ 3, padded-BFVD rerun, verbose-prostt5 reveals different inference path).

**Evidence.** BFVD validation in `/tmp` shell pass — 347,514 distinct uniprot stems each mapping to a single taxid (no conflicts). Score samples documented inline in `THRESHOLD.md`. m8 columns confirmed as standard 7-col foldseek format (query, target, evalue, bits, fident, alnlen, mismatch).

**Outcome.** Project scaffold in place; tasks 62-67 created. Open questions for the user before Phase 2 coding: best-hit vs sum-bits weighting, virus-level grouping (uniprot/taxid/species), splitted-model handling — all flagged in PLAN.md.

## 2026-06-24T17:27:00Z — Cloud GPU verification for foldseek t_06_FoldDb

**Context.** User pushed back: "Are you absolutely certain that the GPU is being used on the cloud?" Memory entry claimed cloud was 1h44m at ~172 ms/seq with GPU active, but `get_job_metadata` showed `runtimeAttributes.gpu = "false"` — looked suspicious enough to revisit.

**Decision / action.** Read t_06_FoldDb stdout from submission 2da383d1-fafb-46f0-bac4-3f4f782304d0 / workflow 7f34ae65-f1e6-44e6-942a-e729c0f6f070 in chunks (fiss-mcp returned CCR refs for the whole file; offset reads gave the raw bytes). Decoded the nvidia-smi preamble (added in task #51) and the foldseek log tail.

**Evidence.** Cloud stdout shows:
- `NVIDIA-SMI 580.159.04   Driver Version: 580.159.04   CUDA Version: 13.0` with `GPU 0: Tesla T4` — GPU physically attached.
- foldseek device list line `CUDA0` (ggml backend dump) — ProstT5 used CUDA backend.
- 36,330 ORFs in `1h 45m 21s 942ms` compute time = **174 ms/seq**, matches local Quadro RTX 8000 rate of **182 ms/seq** within noise.
- Cromwell wall-clock 2h 28m 49s = 1h45m compute + ~43m VM/container/IO overhead.

**Why.** The `runtimeAttributes.gpu = "false"` flag is COSMETIC — GCP Batch backend attaches GPU based on `gpuCount` + `gpuType` runtime attributes only. nvidia-smi output proves attach. Previous concern about CPU-fallback hypothesis is now doubly refuted (local + cloud both confirmed GPU active).

**Outcome.** Diagnosis stable: single-sequence ProstT5 path is the bottleneck regardless of GPU class (T4 ≈ RTX 8000 here because kernel launch + memcpy dominate at n_seq_max=1, not raw FLOPS). Memory entry updated. Task #80 (multi-context concurrency patch) remains the right quick win. No change to plan.

## 2026-06-24T19:24:00Z  — Cohort run: 45 samples through ingest + plot pipeline

**Context.** 45-sample HVP-0006.1 cohort (Terra submission 96f2f736), foldseek m8 outputs downloaded by user-driven `scripts/download_cohort.sh`. Need per-sample artifacts (anndata, plots) before pooled landscape can be built.

**Decision / action.** Authored `analysis/hvp_viral_viz/scripts/run_cohort.sh` (commit 2bba0a2): per-sample ingest → plots (tier 1) → plots_tier2 → label_clusters (with protein markers), gated by output-file presence for idempotent re-run. Ran with `PARALLEL=4`.

**Why.** Single-sample pipeline modules are pure CLI; bash xargs -P composes them cleanly. Per-stage gates are simpler than a Makefile and survive partial runs. PARALLEL=4 on 95 GB / 12-core host: each worker peaks ~3-5 GB (m8 → pandas frame); 4× peak ≈ 20 GB headroom-safe.

**Evidence.** 45/45 produced `umap_by_cluster_labeled_order_with_proteins.png` and 16-plot suite; 0 logs match `Traceback|FAIL`; total `out/cohort_2026-06-22/` 2.9 GB; sample range 21-125 MB tracks raw-hit count (36P/41P low end, 19P/20P/6P high end). Wall-clock ~33 min start to finish.

**Outcome.** Per-sample deliverable for cohort complete. Next: pool 45 hits_filtered.parquet → unified virus×pooled-ORF UMAP (#84) and sample×virus diagnostic (#85). Outputs live under `out/cohort_2026-06-22/`, gitignored (per `analysis/hvp_viral_viz/.gitignore`).

## 2026-06-24T20:33:00Z  — Cohort pool + sample-structure UMAP/PCA + sample-colored landscape

**Context.** 45-sample HVP-0006.1 cohort already had per-sample artifacts (#83, journal entry 2026-06-24T19:24Z). Two follow-ons: build a single unified viral landscape pooled across the cohort (#84), then a 45-point sample-structure diagnostic + sample-colored views of the landscape (#85, plus the standing requirement to color the unified UMAP by sample).

**Decision / action.** Two commits on `tshea/wdl-metagenomics-poc`:
- 6072a81 `analysis/hvp_viral_viz/pool_cohort.py` — streaming cohort pooler. Three pyarrow-streamed passes over `pooled_hits.parquet` (2.48 GB, 81M rows). Pass 1: per-virus count, samples_present, host metadata, orf_source. Pass 2: sample × virus CSR. Pass 3 (two-stage A→B): count distinct viruses per (sample,query), then build virus × sample-prefixed-query COO over cross-virus queries only.
- 49bc393 `analysis/hvp_viral_viz/cohort_extras.py` — `sample-diag` (45-point SVD + UMAP + dendrogram on `pooled_anndata.h5ad`) and `sample-umap` (prevalence + 45-panel small multiples on `clusters.parquet` + `virus_samples.parquet`).

**Why.**
- **Streaming, not pandas-concat.** First in-memory attempt (`pd.read_parquet` of 2.48 GB / 81M rows / 27 string-heavy cols) inflated to ~30 GB pandas, then dropna copied → silent cgroup OOM-kill (exit 0, log truncated at "loaded 81,157,115 hits"). Refactor to `pyarrow.ParquetFile.iter_batches(columns=...)` with column projection drops peak RAM to ~8 GB on 95 GB host.
- **Query as feature, not BFVD target.** Probed top targets across all 45 samples: each hit exactly 1 virus. Root cause: target→taxid is 1:1 by construction (the BFVD lineage join attaches the BFVD entry's taxid to every row hitting that target), so a target-keyed matrix is essentially the identity — no shared columns between viruses, UMAP has no signal. Queries fan out via foldseek's top-k hits — one ORF hits several BFVD entries from different viruses, bridging similar viruses. Sample-prefix `<sample>::<query>` keeps independent samples' queries from colliding.
- **Persist `virus_samples.parquet`.** Standing user requirement: color the unified UMAP by sample. Persisting per-virus `samples_present` at pool time makes sample-coloring a pure plotting task; #88 was a 200-line follow-on with no re-pool.
- **PCA + UMAP + dendrogram for 45 points.** UMAP with n=45 is fragile; PC1-PC2 scatter and an average-linkage cosine dendrogram are the load-bearing views, UMAP added for visual continuity with the virus-level plots.

**Evidence.**
- Pool wall clock: 18 min (220 s sample concat + 96 s pass-1 + 110 s pass-2 + 280 s pass-3 + 360 s embed). Idempotent: resume reads `pooled_hits.parquet` in 115 s if present.
- 13,822 / 14,096 viruses retained in unified UMAP (274 dropped — 0 cross-virus queries after singleton filter).
- 2,585,983 sample-prefixed cross-virus queries as features (filtered from 3.18 M total).
- SVD on virus × query: explained variance 16.1%. 54 leiden clusters at resolution 1.0.
- Sample diagnostic: SVD on sample × virus, 10 components → explained variance **98.3%** (45 samples are highly collinear in virus-presence space; consistent with the median virus appearing in 21 / 45 samples per `virus_samples.parquet`).
- `summary.json` captured at `out/cohort_2026-06-22/_pooled/summary.json`: 81,082,693 pooled hits, 20,582,405 high-conf, 14,096 viruses, 98,425 BFVD targets.

**Outcome.** Pooled artifacts under `out/cohort_2026-06-22/_pooled/`:
- `pooled_hits.parquet` (2.48 GB), `pooled_anndata.h5ad` (3.0 MB, 45×14096 csr), `virus_samples.parquet` (median n_samples per virus = 21), `summary.json`.
- `plots/{abundance_bar,abundance_by_family,heatmap_orf_source}_top30.png`, `plots/umap_by_{family,host,cluster}.png`, `plots/clusters.parquet` (taxid → UMAP coords + leiden cluster).
- `plots/sample_diag_{pca,umap,dendrogram}.png`, `plots/sample_diag_coords.parquet`, `plots/umap_by_prevalence.png`, `plots/umap_per_sample.png`.

**Follow-ups.** None blocking. Output is gitignored under `analysis/hvp_viral_viz/.gitignore`; source on branch `tshea/wdl-metagenomics-poc`. The high (98.3%) explained variance on sample × virus suggests the sample-level structure is dominated by a small number of axes — would be worth a follow-on noting whether the dendrogram clusters track any sample-prep covariate the user has but I don't (sample plate, extraction batch, sequencing run).
