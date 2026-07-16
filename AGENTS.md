# AGENTS.md

**Read this entire file before doing any work in this repo.** Then read [CONTRIBUTING.md](CONTRIBUTING.md) for setup and PR workflow.

Many file edits in this repo look trivial but have non-obvious CI consequences (auto-versioning, manifest format, Trivy filtering, release tagging). This document captures the conventions a fresh agent or contributor would otherwise have to derive by reading two GitHub Actions workflows, two CI shell scripts, and every per-image dotfile.

---

## What this repo is

HVP-LR (Human Virome Project — Long Reads) is a template-style monorepo for the team's long-read bioinformatics work. It holds two kinds of deliverables:

- **Container images** under `docker/<image>/`. One subdirectory per image. Each is independently built, scanned, and published by CI.
- **WDL workflows** under `wdl/*.wdl`. They consume the published images.

Plus the supporting CI plumbing to build, scan, version, publish, and release them.

## Layout

```
HVP-LR/
├── docker/                       # One subdirectory per container image
│   └── <image>/
│       ├── Dockerfile
│       ├── Makefile              # Local build/push (manual); also defines VERSION
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
│   ├── docker.yml                # Build / push / Trivy-scan automation
│   └── cd.yml                    # Release-cutting on main
├── .pre-commit-config.yaml       # WDL validation hook
├── dev-requirements.txt          # Local dev tooling (miniwdl, pre-commit)
├── VERSION                       # Repo-level version bumped by cd.yml on each release
├── LICENSE                       # MIT License (bundled into every image)
├── README.md                     # Human-facing orientation
├── CONTRIBUTING.md               # Setup + PR workflow (humans + agents)
├── CLAUDE.md                     # Pointer for Claude sessions
└── AGENTS.md                     # This file
```

## CI architecture overview

Two workflows, chained:

**`docker.yml` — Docker Build, Push, Scan**

- Triggers: push to **any branch**, PR into `main`, weekly Monday 06:00 UTC cron, manual `workflow_dispatch`. Tag pushes are intentionally NOT a trigger (the release flow creates them; reacting would loop).
- Jobs (push to a branch): `detect-changes` → `bump-versions` (default branch only) → `build-per-arch` (matrix `image × {amd64, arm64}`, each on its native runner — `ubuntu-latest` for amd64, `ubuntu-24.04-arm` for arm64) → `merge-and-scan` (matrix per image: `docker buildx imagetools create` stitches the per-arch intermediate tags into a multi-arch Docker v2 manifest list, then Trivy scans the merged tag).
- Jobs (PR into `main`): `detect-changes` → `pr-build-scan` (amd64-only build + load + Trivy, no push). arm64 is intentionally skipped on PRs — the multi-arch build only runs after merge to `main`.
- Jobs (cron/dispatch): `detect-changes` → `rescan-latest` (pulls and rescans `:latest` for each image, no build).
- **Why native arm runners, not QEMU**: cross-building arm64 under QEMU on an amd64 runner was an order-of-magnitude slower (an HVP-monolith arm64 build alone ran close to the 6-hour job cap). `ubuntu-24.04-arm` runs the arm64 build natively. The trade-off is that each on-push build now publishes a transient `:VERSION-amd64` and `:VERSION-arm64` intermediate tag to GHCR alongside the final multi-arch tags. These intermediates point at the same per-arch digests the manifest list references — they're harmless but visible in the package UI.

**`cd.yml` — CD: Cut release on main**

- Trigger: `workflow_run` listening for `Docker Build, Push, Scan` completing on `main`. Runs only if the upstream conclusion was `success`. Also runs on manual `workflow_dispatch`.
- Bumps repo-root `VERSION`, generates release notes, commits with `[skip ci]`, tags `hvp-lr_v<X.Y.Z>`, creates a GitHub Release.

**Loop prevention.** Auto-bump commits (per-image VERSION and repo-root VERSION) end with `[skip ci]`. The `docker.yml` `detect-changes` job short-circuits on `[skip ci]`. Tag pushes are excluded from `docker.yml` triggers. Both workflows use named `concurrency` groups to serialize runs.

## Docker images — full mechanics

### Per-image required files

Every image directory MUST have all six files. The CI workflow hardcodes paths and will fail loudly if any are missing:

| File | Purpose |
|------|---------|
| `Dockerfile` | Base + install steps. End with `COPY LICENSE /opt/hvp-lr/LICENSE`. |
| `Makefile` | Defines `IMAGE_NAME`, `VERSION`, and inherits build/push targets. Copy verbatim from `docker/hvp-monolith/Makefile` and edit only the two image-specific variables. |
| `env.yaml` | Pinned conda/pip/apt deps (or equivalent for the chosen base image). |
| `.dockerignore` | At minimum: `.git/`, `.github/`, `README.md`, `scripts/`, IDE files, local logs. |
| `.trivyignore` | Per-CVE exceptions with written justification. Start empty. |
| `.trivy-ignore-policy.rego` | Class-based filter (CVSS-vector-based). Today: copy from `docker/hvp-monolith/`. |

The `LICENSE` is staged into the build context at build time (by Makefile and CI), not committed per-image. It's gitignored under `docker/*/LICENSE`.

### Per-image Makefile required variables

Only two variables are image-specific. Everything else is shared boilerplate to copy verbatim.

| Variable | Required | Notes |
|----------|----------|-------|
| `IMAGE_NAME` | **yes** | Repository-relative image name. Must match the directory name under `docker/`. The CI workflow and the auto-bump key off this. |
| `VERSION` | **yes** | Semver `X.Y.Z`. Bumped automatically by CI on `main` if you don't bump it yourself (see Auto-versioning). |
| `REPOS` | inherit | Registry prefixes for `make push`. Standard:<br>`ghcr.io/broadinstitute/hvp-lr`<br>GCR mirror (`us.gcr.io/<GCR_PROJECT>`) is added by CI only when the `GCP_SA_KEY` secret is set — do not include it in the Makefile. |
| `PLATFORMS` | inherit | Default `linux/amd64,linux/arm64`. Override only if a dep genuinely doesn't build on one arch — flag this in PR review. |
| `BUILDX_ATTEST_OFF`, `MANIFEST_OUTPUT` | inherit | Force Docker v2.2 manifest mediatypes. **Do not change** — see Manifest format. |
| `BUILDER` | inherit | Name for the on-demand multi-arch buildx builder. |
| `TAGS`, `TARGETS`, `TAG_ARGS` | inherit | Computed; don't edit. |

Don't touch any of the `.PHONY` recipes (`all`, `build`, `build_no_cache`, `push`, `clean`).

### Tag scheme

Every successful publish of `<image>` carries:

| Tag | When applied |
|-----|--------------|
| `:<VERSION>` | Always. Immutable. |
| `:main` | Push to default branch only. |
| `:latest` | Push to default branch only. |
| `:branch-<name>` | Push to any non-default branch. |
| `:pr-<N>` | Computed for PR builds; NOT pushed (PR builds don't push to a registry). |
| `:<gittag>` | Computed but never triggered — `docker.yml` ignores tag pushes. |

### Auto-versioning rules

`scripts/ci-bump-image-versions.sh` runs only on default-branch pushes. For each changed image:

- If the `VERSION =` line in `docker/<image>/Makefile` was edited in the pushed range → developer set the version, leave it alone.
- If not → patch-bump `X.Y.Z → X.Y.(Z+1)`, commit as `github-actions[bot]` with message `chore(docker): bump versions for <list> [skip ci]`, push.

Practical implications:

- Bump manually if you want `0.0.5 → 0.1.0` or `1.2.3 → 2.0.0`. Otherwise let CI patch-bump.
- Any explicit edit to the `VERSION` line is taken as authoritative — even if it's the same string. **Don't edit `VERSION` cosmetically.**
- Non-default-branch pushes never auto-bump.

### Manifest format

All images MUST be published with the **Docker v2.2 manifest list** mediatype (`application/vnd.docker.distribution.manifest.list.v2+json`), not the OCI image index. Downstream consumers reject OCI mediatypes.

Three settings enforce this in both the per-image Makefile and CI:

- `--provenance=false`
- `--sbom=false`
- `--output type=image,oci-mediatypes=false,push=true`

Provenance and SBOM are off because **either** attestation forces buildx back to OCI mediatypes. If signing/SBOM is needed later, layer it on *after* publishing (e.g. cosign on the v2 manifest).

### Where images publish

| Registry | When | Path |
|----------|------|------|
| GHCR (`ghcr.io`) | Always (non-PR) | `ghcr.io/<owner>/<repo>/<image>:<tag>` (lowercased) |
| GCR (`us.gcr.io`) | Only if `GCP_SA_KEY` secret is set | `us.gcr.io/<GCR_PROJECT>/<image>:<tag>` (project from repo var `GCR_PROJECT`, default `broadinstitute`) |

PRs build but do not push. PR builds are amd64-only and `--load`ed locally as `<image>:scan` for the Trivy step.

GHCR packages are private by default. Visibility is set per-package at `https://github.com/<owner>/<repo>/pkgs/container/<repo>%2F<image>` → Package settings.

### Local build/push targets

```bash
cd docker/hvp-monolith
make build           # Single-arch linux/amd64; loads into local docker as <image>:<VERSION> + :latest
make build_no_cache  # Same as build but --no-cache
make push            # Multi-arch; pushes to every entry in REPOS as Docker v2.2 manifest list
```

`make push` self-bootstraps a `docker-container` builder on first run. It requires `docker login` to each registry in `REPOS`.

`make build` does NOT use the multi-arch builder — it's intentionally single-arch + `--load`ed for fast local iteration.

## Trivy filtering — three independent layers

A vulnerability is suppressed if **any** of the three matches. They are additive, not hierarchical.

1. **Workflow CLI flags** (same for every image, set in `docker.yml`):
   - `severity: CRITICAL,HIGH`
   - `ignore-unfixed: true`
2. **`.trivyignore`** — per-image, per-CVE-ID exceptions with required justification. Format: `CVE-YYYY-NNNNN` per line, `#` for comments.
3. **`.trivy-ignore-policy.rego`** — per-image, class-based rules keyed off the CVSS vector. Today the file is copied verbatim per image; it is the entire class-based policy for that image (there is no shared base policy). See `docker/hvp-monolith/.trivy-ignore-policy.rego` for the canonical content (sections cover AV:P, AV:A, AV:L+UI:R, AV:L+PR:H, AV:L+S:U, availability-only-with-S:U). Quarterly review.

Both `.trivyignore` and `.trivy-ignore-policy.rego` paths are hardcoded in `docker.yml`. **There is no fallback.** Adding an image without these two files makes the Trivy step fail.

Trivy runs twice per image: once as a build gate (table format, `exit-code: 1`), once as SARIF upload to `Security → Code scanning` under category `trivy-<image>`. The weekly cron rescans `:latest` of every image without rebuilding, surfacing newly published CVE data.

**Scan profile is auto-chosen by image size.** `scripts/ci-trivy-profile.sh` inspects the loaded image's `.Size` and picks:

- **Large image** (over `TRIVY_LARGE_IMAGE_BYTES`, default `2147483648` = 2 GiB) — `scanners: vuln`, `timeout: 20m`. Drops the secret scanner, which is too slow on multi-GB scientific containers (it walks every text-ish file and chokes on huge reference DBs / model weights / generated headers) and is intended for source repositories with API keys anyway.
- **Small image** (≤ threshold) — `scanners: vuln,secret`, `timeout: 10m`. Full default scan; small helper images aren't downgraded.

Override the threshold per repo via the `TRIVY_LARGE_IMAGE_BYTES` GitHub Actions repository variable.

## WDL style + validation

Style rules (file layout, boilerplate, task skeleton, naming, output ordering, `t_NN_` call aliasing, etc.) live in [docs/WDL_STYLE_RULES.md](docs/WDL_STYLE_RULES.md). Read it before writing or reformatting any `.wdl` file.

The `miniwdl check --strict` pre-commit hook is the **only** WDL gate. There is no CI safety net. If a contributor commits without running pre-commit (`git commit -n`), broken WDLs land on `main`.

To suppress a known-benign miniwdl lint inline, append `# !LintName` on the offending line, e.g. `# !FileCoercion`.

To run manually on the whole tree:

```bash
pre-commit run --all-files
```

## Release flow

`cd.yml` fires on `workflow_run` after `docker.yml` completes successfully on `main` (and on manual dispatch). It:

1. Reads repo-root `VERSION` (`X.Y.Z`), patch-bumps to `X.Y.(Z+1)`, computes tag `hvp-lr_v<new>`.
2. Builds `release_notes.md` containing:
   - Commit log since the most recent `hvp-lr_v*` tag (or last 50 if none).
   - A table of every `docker/<img>/` showing its current VERSION and the `ghcr.io/...` pull URL.
3. Commits the new `VERSION` as `[skip ci]` and pushes.
4. Creates and pushes annotated git tag `hvp-lr_v<new>`.
5. Creates a GitHub Release (`softprops/action-gh-release`).

Concurrency: `release-main` group, no cancel-in-progress.

## Hard rules — do not

- **Don't bump `VERSION` cosmetically** (per-image or repo-root). Any explicit edit is taken as authoritative.
- **Don't push images by hand from a laptop** outside of testing. CI is the source of truth. A laptop push of `0.0.7` followed by a CI build of `0.0.8` for the same content is confusing for users.
- **Don't change** `--provenance=false`, `--sbom=false`, or `oci-mediatypes=false` on any image's Makefile or in `docker.yml`. Downstream consumers reject OCI mediatypes.
- **Don't skip pre-commit hooks** (`-n` / `--no-verify`). The hook is the only WDL gate.
- **Don't add a new image without all six per-image files**, even if `.trivyignore` and `.trivy-ignore-policy.rego` are empty/copy-pasted. CI hardcodes the paths.
- **Don't edit or amend the bot's `[skip ci]` commits.** They're load-bearing for loop prevention.
- **Don't add CI for WDLs without checking with maintainers first.** There may be a plan; `docs/` is referenced by README but not yet created.
- **Don't disable arm64 multi-arch silently.** If a dep doesn't build on arm64, flag it in PR review and add a comment explaining the opt-out.

## Adding a new image — checklist

1. `mkdir docker/<name>/` (must match `IMAGE_NAME`).
2. Copy these files from `docker/hvp-monolith/` and edit:
   - `Dockerfile` — base + install. End with `COPY LICENSE /opt/hvp-lr/LICENSE`.
   - `Makefile` — change `IMAGE_NAME` and `VERSION = 0.0.1`. Leave everything else.
   - `env.yaml` — pinned deps for the chosen base image (or replace with whatever your base consumes).
   - `.dockerignore` — start from the existing one.
   - `.trivyignore` — keep header; leave entries empty until needed.
   - `.trivy-ignore-policy.rego` — copy verbatim. Edit only if the image violates the assumptions documented at the top of the file (e.g., runs as root, exposes a port).
3. Build locally: `cd docker/<name> && make build` and smoke-test the image.
4. Open a PR. CI runs build + Trivy scan (no push). On merge to `main`, CI auto-bumps if needed and publishes.
5. After first publish, set GHCR package visibility per the team's policy.

## Adding a new WDL — checklist

0. Read [docs/WDL_STYLE_RULES.md](docs/WDL_STYLE_RULES.md) — file layout, task skeleton, naming, call aliasing.
1. Write under `wdl/<workflow>.wdl`.
2. (Optional) Inputs JSON alongside as `wdl/<workflow>.inputs.json`.
3. Run `miniwdl check --strict wdl/<workflow>.wdl` locally. The pre-commit hook runs the same.
4. Reference any required image as `ghcr.io/<owner>/<repo>/<image>:<VERSION>` in the task `runtime { docker: ... }` block. Pin `:VERSION`, not `:latest` or `:main`, for reproducibility.
5. Commit. PR. Merge. (No CI gate — the pre-commit hook is your safety net.)

## Suppressing a Trivy CVE

When a CVE has no upstream fix or no realistic exploit path:

1. Confirm `ignore-unfixed: true` and the rego policy don't already cover it (rerun the Trivy step or run `trivy image` locally).
2. If still surfaced, add to `docker/<image>/.trivyignore`:
   ```
   # CVE-YYYY-NNNNN: <one-line summary> — <why it's accepted: e.g.
   #   "no upstream fix; we never call the affected ctor; review 2026-Q3">
   CVE-YYYY-NNNNN
   ```
3. The `.trivyignore` header notes a quarterly re-evaluation contract.

If the issue is a *class* of CVEs (not a single ID), add a rule to the per-image `.trivy-ignore-policy.rego` instead, with a comment block explaining the rationale (mirror the existing sections' style).

## Pointers

- [README.md](README.md) — short human-facing orientation.
- [CONTRIBUTING.md](CONTRIBUTING.md) — setup + PR workflow shared by humans and agents.
- [LICENSE](LICENSE) — MIT License.
- `.github/workflows/docker.yml` and `.github/workflows/cd.yml` — the source of truth for everything described above. If this file ever drifts from those workflows, the workflows win.
- `scripts/ci-detect-changed-images.sh`, `scripts/ci-bump-image-versions.sh` — CI internals; do not invoke directly.
