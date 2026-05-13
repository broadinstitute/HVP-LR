# HVP-LR

Human Virome Project — Long Reads. WDL pipelines, container images, and supporting scripts for processing PacBio HiFi data and downstream QC.

## Repository layout

```
HVP-LR/
├── docker/                       # One subdirectory per container image
│   └── <image>/                  #   (Dockerfile, Makefile, env.yaml, .dockerignore,
│                                 #    .trivyignore, .trivy-ignore-policy.rego)
├── wdl/                          # WDL workflow definitions
├── scripts/                      # CI helper scripts
├── .github/workflows/            # docker.yml (build/push/scan), cd.yml (release)
├── .pre-commit-config.yaml       # WDL validation hook
├── dev-requirements.txt          # Local dev tooling
├── VERSION                       # Repo-level version
├── LICENSE                       # PolyForm Shield 1.0.0
├── README.md                     # This file
├── CONTRIBUTING.md               # Setup + PR workflow
├── CLAUDE.md                     # Pointer for Claude sessions
└── AGENTS.md                     # CI / Docker / WDL mechanics (read first if contributing)
```

**Where to put new things:**

| What | Where |
|------|-------|
| A new Docker image | `docker/<image-name>/` |
| A new WDL workflow or task | `wdl/<workflow>.wdl` |
| A workflow's inputs JSON | `wdl/<workflow>.inputs.json` (alongside the WDL) |
| Helper / one-off scripts | `scripts/` (avoid if not CI-related; otherwise put inside the relevant `docker/<image>/`) |

## Where to read next

- **Setting up locally / opening a PR** → [CONTRIBUTING.md](CONTRIBUTING.md)
- **CI mechanics, Docker conventions, releases, Trivy filtering, hard rules** → [AGENTS.md](AGENTS.md)
- **License terms** → [LICENSE](LICENSE)

If you're an AI agent (Claude or otherwise), start with [CLAUDE.md](CLAUDE.md) / [AGENTS.md](AGENTS.md).

## Current images

| Image | Subdir | Purpose |
|-------|--------|---------|
| `hvp-monolith` | `docker/hvp-monolith/` | All-in-one QC/alignment toolbox (FastQC, MultiQC, samtools, biopython, pysam, miniwdl, pigz, zstd, GNU parallel, jq) on `mambaorg/micromamba:2.4.0-ubuntu24.04` |

Images are published to `ghcr.io/broadinstitute/hvp-lr/<image>:<version>` and (when GCR is configured) `us.gcr.io/broadinstitute/<image>:<version>`. GHCR packages are private by default.

## License

PolyForm Shield 1.0.0 — see [LICENSE](LICENSE). Source-available, not OSI open source. Any use is permitted except providing a product that competes with HVP-LR or with software the Broad Institute provides using HVP-LR. The LICENSE file is bundled into every built container at `/opt/hvp-lr/LICENSE`.

If you intend to redistribute an image or its derivatives, propagate the LICENSE and any `Required Notice:` lines per the Shield "Notices" section.
