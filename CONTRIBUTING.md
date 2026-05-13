# Contributing to HVP-LR

This document covers setup and the PR workflow for both humans and AI agents working in this repo. For the full mechanics of CI, Docker images, releases, and Trivy filtering, see [AGENTS.md](AGENTS.md). Read AGENTS.md before making non-trivial changes.

## Local dev setup

```bash
python3 -mvenv venv && . venv/bin/activate
pip install -r dev-requirements.txt
pre-commit install
```

The first `pre-commit` run downloads its hooks; subsequent commits are fast.

## Pre-commit hook

A single hook is configured in `.pre-commit-config.yaml`:

- **`miniwdl-check --strict`** — every `.wdl` file is type-checked. The hook refuses commits that introduce strict-mode WDL warnings or errors.

This hook is **the only WDL gate**. There is no CI safety net. Don't bypass it (`-n` / `--no-verify`).

To suppress a known-benign miniwdl lint, append `# !LintName` on the offending line, e.g. `# !FileCoercion`.

To run the hook on the whole tree:

```bash
pre-commit run --all-files
```

## Branching and PRs

- Branch from `main`. Use a descriptive prefix: `feat/...`, `fix/...`, `docs/...`, `ci/...`.
- Open the PR against `main`.
- CI runs build + Trivy scan for every changed Docker image (no push for PRs). Build failures or Trivy `HIGH/CRITICAL` findings block merge.
- After merge to `main`, CI auto-bumps changed images, publishes to GHCR (and GCR if configured), and `cd.yml` cuts a release tag.

## Commit messages

- Be descriptive. The body matters more than the subject for non-trivial changes.
- **Avoid the `chore(docker):` and `chore(release):` prefixes** — those are reserved for the `github-actions[bot]` `[skip ci]` commits. Using them in human commits makes the history harder to scan.
- **Don't include `[skip ci]` in human commits.** It's load-bearing for CI loop prevention; using it elsewhere will silently skip CI on your changes.

## Common contribution paths

### Adding a Docker image

1. Create `docker/<name>/` with all six required files (Dockerfile, Makefile, env.yaml, .dockerignore, .trivyignore, .trivy-ignore-policy.rego).
2. Build and smoke-test locally: `cd docker/<name> && make build`.
3. PR. Merge. CI handles versioning and publishing.

The CI workflow hardcodes paths — missing files will fail the build. See [AGENTS.md → Adding a new image](AGENTS.md#adding-a-new-image--checklist) for the full checklist and the reasoning behind each file.

### Adding a WDL workflow

1. Write under `wdl/<workflow>.wdl`. Pin docker image references to `:<VERSION>`, not `:latest` or `:main`.
2. Run `miniwdl check --strict wdl/<workflow>.wdl` locally (the pre-commit hook runs the same).
3. PR. Merge.

### Suppressing a Trivy CVE

When a CVE has no upstream fix or no realistic exploit path against this image:

1. Add the CVE ID to `docker/<image>/.trivyignore` with a one-line justification:
   ```
   # CVE-YYYY-NNNNN: <one-line summary> — <why accepted; review date>
   CVE-YYYY-NNNNN
   ```
2. Quarterly re-evaluation is part of the contract. Remove the entry when the underlying package ships a fix.

For *classes* of CVEs (not single IDs), see [AGENTS.md → Trivy filtering](AGENTS.md#trivy-filtering--three-independent-layers) for the per-image `.trivy-ignore-policy.rego`.

## Issues and questions

File issues at https://github.com/broadinstitute/HVP-LR/issues. For sensitive security issues, do not file a public issue — contact the maintainers directly.

## License

Contributions are accepted under [PolyForm Shield 1.0.0](LICENSE) (source-available, not OSI open source). The license file is bundled into every built container at `/opt/hvp-lr/LICENSE`.
