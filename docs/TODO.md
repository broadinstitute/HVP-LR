# TODO / backlog

Durable notes for work we've deliberately deferred. Add items here rather than
losing them in commit messages or chat.

## CI: run tests (WDL + Python)

**There is no CI test job.** `docker.yml` builds + Trivy-scans images; `cd.yml`
tags/releases. Neither runs any test. WDL correctness rests solely on the
`miniwdl check --strict` pre-commit hook (see AGENTS.md §"only WDL gate"), which
a contributor can bypass with `git commit -n`, and which only lints — it does not
*run* anything.

Tests that exist but nothing enforces:
- `tests/python/` — unit tests mirroring WDL task Python heredocs.
- `docker/x_dosage/tests/` — 18 pytest for the x_dosage classifier/calibration.
- `tests/wdl/test_sexchrom_karyotype.py` — end-to-end SexChromKaryotype run
  (self-skips without docker/miniwdl/image).

**Add a CI workflow that runs these on PRs.** Suggested, cheapest-first:
1. `miniwdl check --strict` on all `wdl/**/*.wdl` (mirror the pre-commit hook so
   `-n` bypass can't land broken WDL on main). Cheap, no docker.
2. `pytest docker/x_dosage/tests tests/python` on relevant path changes. Cheap,
   no docker. Would have caught the x_dosage config-path breakage automatically.
3. (Optional, heavier) the `tests/wdl/` integration test on a docker-enabled
   runner with the built image — real pipeline execution.

Owner: unassigned. Raised 2026-08 during x_dosage WDL/image work.
