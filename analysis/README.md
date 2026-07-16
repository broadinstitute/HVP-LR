# analysis/

Tools and analyses that operate on data produced by the HVP-LR pipelines.

This folder is **not** part of the WDL workflow execution path. WDL tasks
live under [`wdl/`](../wdl/) and container images under [`docker/`](../docker/).
What lives here instead is downstream code that consumes pipeline outputs
(hits TSVs, AnnData files, BFVD reference joins, foldseek m8 tables, etc.)
and turns them into figures, reports, cluster labels, marker tables, and
other interpretive artifacts.

## When to add something here

Put a new project under `analysis/<name>/` when **all** of these are true:

- It reads, transforms, plots, or annotates data emitted by an HVP-LR
  workflow (or a sibling reference dataset like BFVD / ICTV).
- It does **not** need to run inside a Cromwell task as part of a sample's
  primary pipeline. If it needs to run on every sample as part of the
  pipeline, it belongs in `wdl/` with a container in `docker/`.
- It is meaningful at the project level — useful to more than one sample,
  more than one user, or more than one experiment. One-off scratch
  notebooks do not belong here; check them into your own workspace.

## Conventions

Each project lives in its own subdirectory and is self-contained:

```
analysis/<name>/
├── README.md             # what it does, how to run it, inputs/outputs
├── pyproject.toml        # if Python, for `pip install -e .`
├── <package>/            # source (or src/<package>/ — either layout is fine)
├── tests/                # pytest-style unit tests
└── docs/                 # any longer-form design docs
```

Rules:

- **No bulk data in git.** Reference parquets, AnnData files, foldseek m8
  outputs, per-sample plots, and UniProt caches do **not** live here. Stash
  them in the relevant GCS bucket and document the URI in the project's
  README. Add `*.parquet`, `*.h5ad`, `out/`, `refs/cache/`, etc. to a local
  `.gitignore` inside the project.
- **Pin your dependencies.** A `pyproject.toml` (or `requirements.txt` with
  hashes) keeps the analysis reproducible.
- **Document the upstream contract.** State which workflow's outputs you
  consume and which columns / schema you assume. If the upstream schema
  changes, the analysis must update in the same PR.
- **Tests are not optional.** At minimum a smoke test per module that
  loads a tiny fixture and runs the main entry point.

## Current projects

| Project | Path | Purpose |
|---------|------|---------|
| `hvp_viral_viz` | [`analysis/hvp_viral_viz/`](hvp_viral_viz/) | Downstream visualization of `HvpViralProteinAnnotation` outputs: virus×ORF embedding, leiden clustering with rank-coherent resolution selection, biologically-labelled UMAPs, per-cluster marker proteins from UniProt. Consumes `hits.tsv` + BFVD reference joins. |
