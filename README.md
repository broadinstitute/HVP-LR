# HVP-LR

Human Virome Project — Long Reads. WDL pipelines, container images, and supporting scripts for processing PacBio HiFi data and downstream QC.

## Repository layout

```
HVP-LR/
├── docker/                       # One subdirectory per container image
│   └── <image>/
│       ├── Dockerfile
│       ├── Makefile              # Local build/push (manual)
│       ├── env.yaml              # micromamba/conda environment (or equivalent)
│       ├── .dockerignore
│       ├── .trivyignore          # Per-CVE exceptions with justifications
│       └── .trivy-ignore-policy.rego  # Rego policy for class-based filtering
├── wdl/                          # WDL workflow definitions
│   └── *.wdl
├── scripts/                      # CI helper scripts (do not invoke directly)
│   ├── ci-detect-changed-images.sh
│   └── ci-bump-image-versions.sh
├── .github/workflows/
│   └── docker.yml                # Build / push / Trivy-scan automation
├── docs/                         # Internal docs (e.g. CI testing plan)
├── .pre-commit-config.yaml       # WDL validation hook
├── dev-requirements.txt          # Local dev tooling (miniwdl, pre-commit)
├── LICENSE                       # PolyForm Shield 1.0.0
└── README.md
```

**Where to put new things:**

| What | Where |
|------|-------|
| A new Docker image | `docker/<image-name>/` |
| A new WDL workflow or task | `wdl/<workflow>.wdl` |
| A workflow's inputs JSON | `wdl/<workflow>.inputs.json` (alongside the WDL) |
| Helper / one-off scripts | `scripts/` (avoid if not CI-related; otherwise put inside the relevant `docker/<image>/`) |

## Dev environment setup

```bash
python3 -mvenv venv && . venv/bin/activate && pip install -r dev-requirements.txt
pre-commit install
```

The first `pre-commit` run downloads its hooks; subsequent commits are fast.

## Pre-commit hooks

Defined in `.pre-commit-config.yaml`. Currently:

- **`miniwdl-check --strict`** — every `.wdl` file is type-checked. The hook will refuse a commit that introduces a strict-mode WDL warning or error. To silence a known-benign miniwdl lint inline, append `# !<LintName>` on the offending line (e.g. `# !FileCoercion`).

To run the hook manually on the whole tree:

```bash
pre-commit run --all-files
```

## Docker images

### Current images

| Image | Subdir | Purpose |
|-------|--------|---------|
| `hvp-monolith` | `docker/hvp-monolith/` | All-in-one QC/alignment toolbox (FastQC, MultiQC, samtools, biopython, pysam, miniwdl, pigz, zstd, GNU parallel, jq) on `mambaorg/micromamba:2.4.0-ubuntu24.04` |

### Manual build & push (developer workflow)

Each image has its own `Makefile` with three useful targets:

| Target | What it does |
|--------|--------------|
| `make build` | Single-arch `linux/amd64` build, `--load`ed into the local docker daemon as `<image>:<VERSION>` and `<image>:latest`. Fast — no cross-arch QEMU. |
| `make build_no_cache` | Same, but `--no-cache`. |
| `make push` | Multi-arch (`linux/amd64,linux/arm64`) buildx build, pushed to every entry in `REPOS` as a **Docker v2.2 manifest list** (`oci-mediatypes=false`). Self-bootstraps a `docker-container` builder on first run. Requires `docker login` to each registry in `REPOS`. |

Local example:

```bash
cd docker/hvp-monolith
make build
docker run --rm -it hvp-monolith:0.0.1 fastqc --version
```

### Adding a new image

1. `mkdir docker/<name>/`
2. Drop in:
   - `Dockerfile` — base + install steps. End with `COPY LICENSE /opt/hvp-lr/LICENSE` (the CI / Makefile stage `LICENSE` into the build context for you).
   - `Makefile` — copy from `docker/hvp-monolith/Makefile` and edit the two image-specific variables (see below).
   - `env.yaml` (or whatever your base image consumes) — pinned conda/pip/apt deps.
   - `.dockerignore` — at minimum, ignore `.git/`, `.github/`, `README.md`, `scripts/`, IDE files, local logs.
   - `.trivyignore` — start empty; add CVE IDs with a written justification only when unavoidable.
   - `.trivy-ignore-policy.rego` — copy from `docker/hvp-monolith/`; tweaks Trivy filtering rules.

3. Commit. The next push that touches `docker/<name>/` triggers CI for the new image automatically (matrix-based, no workflow edits required).

### Required variables in a per-image Makefile

Only two variables are image-specific. Everything else is shared boilerplate that should be copied verbatim from `docker/hvp-monolith/Makefile`.

| Variable | Required | Notes |
|----------|----------|-------|
| `IMAGE_NAME` | **yes** | Repository-relative image name. Must match the directory name under `docker/`. The CI workflow and automatic version-bump key off this. |
| `VERSION` | **yes** | Semver `X.Y.Z`. Bumped automatically by CI on `main` if you don't bump it yourself (see [Auto-versioning](#auto-versioning)). |
| `REPOS` | inherit | List of registry prefixes to push to. Standard value:<br>`us.gcr.io/broadinstitute`<br>`ghcr.io/broadinstitute`<br>Adjust only if this image needs to live somewhere else. |
| `PLATFORMS` | inherit | Default `linux/amd64,linux/arm64`. Override only if a dependency genuinely doesn't build on one arch — flag this in a comment if you do. |
| `BUILDX_ATTEST_OFF`, `MANIFEST_OUTPUT` | inherit | Force Docker v2.2 manifest mediatypes. **Do not change** — see [Manifest format](#manifest-format). |
| `BUILDER` | inherit | Name for the on-demand multi-arch buildx builder. |
| `TAGS`, `TARGETS`, `TAG_ARGS` | inherit | Computed; don't edit. |

You should not need to touch any of the `.PHONY` recipes (`all`, `build`, `build_no_cache`, `push`, `clean`).

### Manifest format

All images in this repo MUST be published with the **Docker v2.2 manifest list** mediatype (`application/vnd.docker.distribution.manifest.list.v2+json`) rather than the OCI image index. Downstream consumers of HVP-LR images reject OCI mediatypes.

The Makefile and CI workflow enforce this via three settings:

- `--provenance=false`
- `--sbom=false`
- `--output type=image,oci-mediatypes=false,push=true`

Provenance and SBOM attestations are off because **both** force buildx to switch back to OCI mediatypes. If you ever want signing/SBOM, you'll need to layer it on after publishing (e.g. cosign on the v2 manifest).

## Continuous integration

Workflow: `.github/workflows/docker.yml`. Runs on every push (any branch + tags), every PR to `main`, weekly on Monday at 06:00 UTC, and on manual dispatch.

### What CI does

1. **Detects which images changed** by diffing `docker/*/` between the push's base and HEAD (`scripts/ci-detect-changed-images.sh`).
2. **Bumps the patch `VERSION`** on `main` if a dev modified an image without bumping it themselves (`scripts/ci-bump-image-versions.sh`). Commits `chore(docker): bump versions for [...] [skip ci]` and pushes. The `[skip ci]` marker prevents the bump commit from retriggering the workflow.
3. **Builds each changed image** with buildx for `linux/amd64,linux/arm64` and pushes to **GHCR** at `ghcr.io/<owner>/<repo>/<image>:<tag>`. If a `GCP_SA_KEY` secret is configured, it also mirrors to `us.gcr.io/<GCR_PROJECT>/<image>:<tag>`.
4. **Scans** the resulting image with Trivy (`CRITICAL,HIGH`, `ignore-unfixed`, per-image `.trivyignore` + `.trivy-ignore-policy.rego`). Findings appear in the repo's **Security → Code scanning** tab. Build fails on unsuppressed HIGH/CRITICAL findings.
5. **Weekly cron / manual dispatch** rescans `:latest` of every image without building, picking up newly published CVE data.

### Tag scheme

For every successful publish of `<image>`:

| Tag | When applied |
|-----|--------------|
| `:<VERSION>` | Always. Immutable. |
| `:main` | Push to default branch only. |
| `:latest` | Push to default branch only. |
| `:branch-<name>` | Push to any non-default branch. |
| `:pr-<N>` | PR builds (PR check only; not currently pushed to the registry). |
| `:<gittag>` | Push of an annotated/lightweight git tag. |

### Auto-versioning

The rule (see `scripts/ci-bump-image-versions.sh` for the exact diff check):

- Dev modifies `docker/<image>/` files **without** changing the `VERSION` line in the Makefile → CI bumps `VERSION` from `X.Y.Z` to `X.Y.(Z+1)`, commits with `[skip ci]`, pushes.
- Dev modifies the `VERSION` line themselves (any bump kind) → CI leaves it alone and publishes the version the dev set.
- Non-default branch pushes never auto-bump; only the default branch.

In practice: bump manually if you want `0.0.5 → 0.1.0` or `1.2.3 → 2.0.0`. Otherwise let CI patch-bump.

## License

PolyForm Shield 1.0.0 — see `LICENSE`. Source-available, not OSI open source. Any use is permitted except providing a product that competes with HVP-LR or with software the Broad Institute provides using HVP-LR. The LICENSE file is bundled into every built container at `/opt/hvp-lr/LICENSE`.

If you intend to redistribute an image or its derivatives, propagate the LICENSE and any `Required Notice:` lines per the Shield "Notices" section.

## Useful pointers for new devs

- **Validate a WDL locally before committing:** `miniwdl check --strict wdl/<file>.wdl`. The pre-commit hook runs the same command.
- **Suppress a known-benign WDL lint** (e.g. unavoidable `String → File` coercion at a scatter source): append `# !FileCoercion` on the line.
- **Don't push images by hand from your laptop unless you're testing** — CI is the source of truth. If your laptop push lands a `0.0.7` and CI later builds `0.0.8` for the same content, that's confusing for users.
- **Don't bump VERSION cosmetically** — CI's diff-aware bump means any explicit edit to the VERSION line is taken as authoritative.
- **When adding a new container, also add per-image `.trivyignore` and `.trivy-ignore-policy.rego`** even if both are empty/copy-pasted. The workflow expects them at fixed paths.
- **arm64 builds use QEMU emulation in CI** — they're slow (often 5–10× longer than amd64). If a dep simply doesn't build on arm64 we'll need to opt that image out of multi-arch; flag it in PR review rather than silently disabling.
- **GHCR packages are private by default.** Visibility is set per-package at `https://github.com/<owner>/<repo>/pkgs/container/<repo>%2F<image>` → Package settings.
- **Trivy failures with no clear fix** — add the CVE ID to `docker/<image>/.trivyignore` with a written justification (the file has a template at the top). Quarterly review of accepted CVEs is part of the contract.
