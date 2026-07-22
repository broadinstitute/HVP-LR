# hvp_viral_viz

Downstream visualization of `HvpViralProteinAnnotation` outputs. Ingests
foldseek m8 hits against BFVD and produces single-cell-style virome
visualizations: virus×ORF count matrices, leiden clustering with
rank-coherent ICTV resolution selection, UMAPs colored by taxonomy / host,
per-cluster marker proteins from UniProt.

This package lives under [`analysis/`](../README.md) — it is **not**
executed by the WDL pipeline. Run it locally (or in a notebook) against
pipeline outputs you have already downloaded from a Terra workspace.

## Layout

```
hvp_viral_viz/
├── README.md         ← this file
├── PLAN.md           ← implementation plan, sequenced
├── THRESHOLD.md      ← bit-score / e-value cutoff + rationale
├── HOST_FILL.md      ← host species fill strategy + alternatives
├── __init__.py
└── *.py              ← package modules (flat layout, package name = folder name)
```

Pending follow-ups: `pyproject.toml` for `pip install -e .`, `tests/` with
smoke tests per module.

Data lives **outside** the repo:

```
<workdir>/refs/          ← processed BFVD + UniProt reference tables (gitignored)
<workdir>/refs/cache/    ← UniProt fetch caches (gitignored)
<workdir>/out/<sample>/  ← per-sample AnnData + plots (gitignored)
```

## Install

No setup.py / pyproject yet. Run modules in place:

```bash
cd analysis
PYTHONPATH=. PYTHONPATH=.. python -m hvp_viral_viz.<module> [args]
```

Heavy deps: `scanpy`, `leidenalg`, `igraph`, `umap-learn`, `anndata`,
`scikit-learn`, `matplotlib`. No GPU required. Install with `pip`:

```bash
pip install anndata scanpy leidenalg igraph umap-learn scikit-learn \
            scipy numpy pandas pyarrow matplotlib requests
```

## Inputs

Per sample, from a `HvpViralProteinAnnotation` run:
- `<sample>.vs_bfvd.m8` — raw foldseek hits (`t_08_Format` output)
- `<sample>.vs_bfvd.annotated_hits.tsv` — annotated best hits (`t_09`)

BFVD reference (one-time):
- `bfvd_metadata.tsv` — uniprot→(model_id, avg_plddt, ptm, splitted, version)
- `bfvd_taxid.tsv` — model→taxid
- `bfvd_taxid_rank_scientificname_lineage.tsv` — model→taxid + GTDB-style lineage

External fill (host species, see [`HOST_FILL.md`](HOST_FILL.md)):
- ICTV Virus Metadata Resource (VMR) — family→host_group
- UniProt REST `organism_host` — per uniprot
- UniProt REST `protein_name` — per uniprot (for cluster-marker labels)
- NCBI Virus host metadata — per taxid

## Pipeline

Numbered entry points live in [`scripts/`](scripts/) and run in order
from the package root. Steps 01–06 build reference tables (one-time per
BFVD release). Steps 07–13 prepare and analyze a cohort.

```bash
cd analysis/hvp_viral_viz

# --- one-time reference build ---
bash scripts/01_fetch_refs.sh --out-dir refs/downloads
bash scripts/02_build_bfvd_refs.sh
bash scripts/03_uniprot_host.sh --source bfvd
bash scripts/04_host_taxid_lookup.sh
bash scripts/05_build_host_table.sh
bash scripts/06_uniprot_protein_name.sh --source bfvd   # opt-in; req'd for label_clusters --with-protein-markers

# --- per-cohort run ---
python scripts/07_terra_json_to_tsv.py --dir data/terra_tables/<date>/    # specimen labels
python scripts/08_build_manifest.py                                       # verify m8 GCS paths
bash   scripts/09_download_cohort.sh                                      # pull m8 files
PARALLEL=4 bash scripts/10_run_cohort.sh                                  # per-sample ingest+plots+labels
bash scripts/11_pool_cohort.sh --cohort-dir out/cohort_<date>
bash scripts/12_cohort_extras.sh --pooled-dir out/cohort_<date>/_pooled all
bash scripts/13_cohort_specimen.sh \
     --pooled-dir out/cohort_<date>/_pooled \
     --table-tsv  data/terra_tables/<date>/HVP-0006_1.tsv
```

Step 10 is idempotent on per-stage marker files; delete a sample's `out/`
subdir to force rebuild. The `--with-protein-markers` flag is opt-in
inside step 10: default `label_clusters.py` output is byte-stable
against runs that pre-date the protein-marker overlay. See
[`PLAN.md`](PLAN.md) for the implementation sequence.

Optional/diagnostic modules (no numbered wrapper):

```bash
PYTHONPATH=.. python -m hvp_viral_viz.scan_resolution --sample-dir out/<sample>
```

## Approach

Single-cell analogy:
- "cells" = viral taxids that survived the hit-quality filter
- "features" = query ORFs (unique `query` IDs)
- "expression" = hit count per (virus, ORF)
- clustering = leiden on a cosine kNN graph over normalized SVD
- labels = dominant ICTV rank + dominant host_group, purity-tagged
- markers = `sc.tl.rank_genes_groups` over the virus×ORF matrix, top hits
  joined to UniProt `protein_name`

Resolution is picked by maximizing ARI of the leiden partition against
ICTV ranks (class / order / family / genus). The resolution+rank with the
highest ARI (NMI as tiebreaker) defines the rank-coherent partition —
clusters at that resolution sit on a single hierarchy level rather than
mixing ranks.

## Inputs the package assumes are already cached

| Path | Producer | Notes |
|------|----------|-------|
| `refs/uniprot_taxid.parquet` | `build_bfvd_refs` | BFVD uniprot→taxid join |
| `refs/taxid_lineage.parquet` | `build_bfvd_refs` | NCBI lineage per taxid |
| `refs/taxid_host.parquet` | `build_host_table` | host_group per virus taxid |
| `refs/cache/uniprot_host.parquet` | `uniprot_host` | UniProt `virus_hosts` field |
| `refs/cache/uniprot_protein_name.parquet` | `uniprot_protein_name` | UniProt `protein_name` field (opt-in) |

The default `--refs-dir` is `/workspace/viral-viz/refs` for historical
reasons (this package was developed standalone); override with the CLI
flag once you settle on a project-relative reference root.

## Status

Functional end-to-end on `HVP-0006.1_34P`. Pending: tests, project-relative
`--refs-dir` default, cohort UMAP once ≥ 2 samples are ingested.
