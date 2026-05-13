# CLAUDE.md

**Before doing any work in this repo, read [AGENTS.md](AGENTS.md) in its entirety.** Then read [CONTRIBUTING.md](CONTRIBUTING.md).

This is not optional and not skippable for "small" or "obvious" changes. Many file edits in this repo look trivial but have non-obvious CI consequences:

- Editing a `Makefile` `VERSION` line — even setting it to the same value — is taken as authoritative by the auto-bump script and overrides CI's patch-bump.
- Adding a new `docker/<image>/` directory without all six required files (Dockerfile, Makefile, env.yaml, .dockerignore, .trivyignore, .trivy-ignore-policy.rego) will fail CI; the workflow hardcodes paths with no fallback.
- Changing `--provenance`, `--sbom`, or `oci-mediatypes` flags will publish images that downstream consumers reject.
- Skipping the `pre-commit` hook (`-n` / `--no-verify`) lets broken WDLs land on `main`; there is no CI safety net.
- Including `[skip ci]` in a human commit message will silently skip CI on that commit.

AGENTS.md documents these and many other conventions. Read it first. If a user request and AGENTS.md disagree, surface the conflict before acting.
