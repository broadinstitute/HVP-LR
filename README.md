# HVP-LR

Human Virome Project — Long Reads. WDL pipelines, container images, and supporting scripts for processing PacBio HiFi data and downstream QC.

## Repository layout

```
HVP-LR/
├── docker/                       # One subdirectory per container image
│   └── <image>/                  #   (Dockerfile, Makefile, env.yaml, .dockerignore,
│                                 #    .trivyignore, .trivy-ignore-policy.rego)
├── wdl/                          # WDL workflow definitions
│   ├── structs/Structs.wdl       #   Shared RuntimeAttr struct
│   ├── tasks/<Category>/         #   Task modules (Preprocessing, QC, Utility, ...)
│   └── pipelines/<Tech>/<Cat>/   #   End-to-end workflows (<Tech> = PacBio | ONT | ILMN | TechAgnostic)
├── analysis/                     # Downstream tools and analyses consuming
│   └── <project>/                #   pipeline outputs (figures, cluster labels,
│                                 #   marker tables, …). NOT executed by WDL.
├── docs/                         # Project docs (WDL style rules, etc.)
├── scripts/                      # CI helper scripts
├── .github/workflows/            # docker.yml (build/push/scan), cd.yml (release)
├── .pre-commit-config.yaml       # WDL validation hook
├── .dockstore.yml                # Dockstore workflow registration manifest
├── dev-requirements.txt          # Local dev tooling
├── VERSION                       # Repo-level version
├── LICENSE                       # MIT License
├── README.md                     # This file
├── CONTRIBUTING.md               # Setup + PR workflow
├── CLAUDE.md                     # Pointer for Claude sessions
└── AGENTS.md                     # CI / Docker / WDL mechanics (read first if contributing)
```

**Where to put new things:**

| What | Where |
|------|-------|
| A new Docker image | `docker/<image-name>/` |
| A new task module | `wdl/tasks/<Category>/<Name>.wdl` |
| A new pipeline / workflow | `wdl/pipelines/<Tech>/<Category>/<Name>.wdl` |
| A workflow's inputs JSON | Alongside the WDL as `<Name>.inputs.json` |
| Helper / one-off scripts | `scripts/` (avoid if not CI-related; otherwise put inside the relevant `docker/<image>/`) |
| A downstream analysis / visualization tool | `analysis/<project>/` — see [analysis/README.md](analysis/README.md) for conventions |

## The `analysis/` folder

[`analysis/`](analysis/) holds downstream tools and analyses that operate on
data produced by HVP-LR pipelines (and sibling reference datasets like
BFVD / ICTV). It is **not** part of the WDL execution path — code here is
invoked by humans (or notebooks) against pipeline outputs to produce
figures, cluster labels, marker tables, and other interpretive artifacts.

Each project lives in its own self-contained subdirectory with its own
`README.md`, dependency manifest, source, and tests. See
[`analysis/README.md`](analysis/README.md) for full conventions and the
list of current projects.

## Where to read next

- **Setting up locally / opening a PR** → [CONTRIBUTING.md](CONTRIBUTING.md)
- **CI mechanics, Docker conventions, releases, Trivy filtering, hard rules** → [AGENTS.md](AGENTS.md)
- **WDL style rules** (file layout, task skeleton, naming, call aliasing) → [docs/WDL_STYLE_RULES.md](docs/WDL_STYLE_RULES.md)
- **License terms** → [LICENSE](LICENSE)

If you're an AI agent (Claude or otherwise), start with [CLAUDE.md](CLAUDE.md) / [AGENTS.md](AGENTS.md).

## Current images

| Image | Subdir | Purpose |
|-------|--------|---------|
| `hvp-monolith` | `docker/hvp-monolith/` | All-in-one QC/alignment toolbox (FastQC, MultiQC, samtools, biopython, pysam, miniwdl, pigz, zstd, GNU parallel, jq) on `mambaorg/micromamba:2.4.0-ubuntu24.04` |

Images are published to `ghcr.io/broadinstitute/hvp-lr/<image>:<version>` and (when GAR is configured) `us-central1-docker.pkg.dev/<GCR_PROJECT>/<GAR_REPO>/<image>:<version>`. GHCR packages are private by default.

## Current workflows

| Workflow | Path | Purpose |
|----------|------|---------|
| `HifiReadQCPipeline` | `wdl/pipelines/PacBio/QC/HifiReadQCPipeline.wdl` | Single-sample PacBio HiFi read QC: taxonomy lookup, BAM→FASTQ, seqkit + kraken2 stats, per-sample human-readable report |

Workflows are registered with Dockstore via [.dockstore.yml](.dockstore.yml).

## License

MIT — see [LICENSE](LICENSE). The LICENSE file is bundled into every built container at `/opt/hvp-lr/LICENSE`. If you redistribute the source or an image, MIT requires the copyright notice and license text to travel with it.
