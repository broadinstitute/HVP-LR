# Testing the Docker CI pipeline

Test sequence to validate `.github/workflows/docker.yml` end-to-end after the
repo is first pushed to GitHub. Run from the repo root unless noted.

## 0. Push repo + enable workflow write perms

```bash
git push -u origin main
```

Then in the GitHub web UI:

**Settings → Actions → General → Workflow permissions** = `Read and write permissions`. Save.

This is required so the `bump-versions` job can commit the auto-bump back to
`main`, and so `build-scan` can push packages to GHCR.

---

## 1. Test A — PR build (no push)

```bash
git checkout -b test/ci-1
printf '\n# noop\n' >> docker/hvp-monolith/Dockerfile
git -c user.email=11667487+jonn-smith@users.noreply.github.com \
    -c user.name="Jonn Smith" \
    commit -am "test: PR build"
git push -u origin test/ci-1
gh pr create --fill --base main
```

**Expect:**
- `detect-changes` outputs `images=["hvp-monolith"]`
- `bump-versions` is skipped (not on main)
- `build-scan` builds a local single-arch image + runs Trivy
- No GHCR push happens
- PR check is green if no HIGH/CRITICAL CVEs

---

## 2. Test B — merge to main → auto-bump + multi-arch push

Merge the PR. Watch the Actions tab.

**Expect:**
- `bump-versions` commits `chore(docker): bump versions for ["hvp-monolith"] [skip ci]` to main
- `build-scan` builds multi-arch and pushes:
  - `ghcr.io/<owner>/hvp-lr/hvp-monolith:0.0.2`
  - `ghcr.io/<owner>/hvp-lr/hvp-monolith:latest`
  - `ghcr.io/<owner>/hvp-lr/hvp-monolith:main`
- The bump commit does **NOT** trigger a second workflow run (verify in Actions tab)
- `Settings → Packages` shows the `hvp-monolith` package

---

## 3. Verify manifest + arches + LICENSE in image

```bash
docker manifest inspect ghcr.io/<owner>/hvp-lr/hvp-monolith:0.0.2 | head -30
docker run --rm ghcr.io/<owner>/hvp-lr/hvp-monolith:0.0.2 cat /opt/hvp-lr/LICENSE | head -5
```

**Expect:**
- `mediaType: application/vnd.docker.distribution.manifest.list.v2+json`
- Two `manifests` entries: `linux/amd64`, `linux/arm64`
- LICENSE output starts with `# PolyForm Shield License 1.0.0`

---

## 4. Test C — non-docker change ignored

```bash
git checkout main && git pull
echo "noop" >> README.md
git -c user.email=11667487+jonn-smith@users.noreply.github.com \
    -c user.name="Jonn Smith" \
    commit -am "test: docs only"
git push
```

**Expect:** `detect-changes` outputs `any=false`; all downstream jobs are skipped.

---

## 5. Test D — manual dispatch rescan

Actions tab → **Docker Build, Push, Scan** → **Run workflow** → `main`.

**Expect:**
- `rescan-latest` pulls `:latest` of every image
- Trivy scans it and uploads SARIF under category `trivy-hvp-monolith`
- No build, no push happens

---

## 6. Test E — dev-bumped VERSION respected

```bash
git checkout -b test/ci-3
sed -i 's/^VERSION = .*/VERSION = 0.1.0/' docker/hvp-monolith/Makefile
printf '\n# noop2\n' >> docker/hvp-monolith/Dockerfile
git -c user.email=11667487+jonn-smith@users.noreply.github.com \
    -c user.name="Jonn Smith" \
    commit -am "test: dev manual bump"
git push -u origin test/ci-3
gh pr create --fill --base main
gh pr merge --merge
```

**Expect:** `bump-versions` runs but reports `committed=false` (the script saw
the VERSION line already changed in the diff). `build-scan` pushes `:0.1.0`,
**not** `:0.1.1`.

---

## Failure modes to watch

- **403 on package push** → workflow write perms not enabled (see Step 0)
- **Bump job cannot push** → same root cause
- **Trivy fails on base-image CVEs** → add CVE IDs to
  `docker/hvp-monolith/.trivyignore` with a written justification
- **More than 1 workflow run per push** → grep the bump commit message for the
  literal substring `[skip ci]`; if it's missing, the recursion guard is broken
