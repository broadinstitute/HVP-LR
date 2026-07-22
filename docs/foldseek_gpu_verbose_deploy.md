# Deploying the verbose foldseek-gpu image

This document covers rolling out
`docker/foldseek-gpu/patches/verbose-prostt5.patch` and the matching
`Dockerfile` wiring to a fresh image, then routing the two HVP-LR
tasks that hit ProstT5 (t_06_FoldDb, t_07_Search) at the new tag.

Companion to `docs/foldseek_gpu_test_plan.md`. Read
`docker/foldseek-gpu/patches/README.md` for what the patch changes.

## TL;DR

1. Open PR with the patch + Dockerfile changes (already committed
   locally on the working branch).
2. Merge to `main`. CI auto-bumps `docker/foldseek-gpu/Makefile`
   VERSION (10.0.1 → 10.0.2) and publishes the multi-arch tag to
   both GHCR and GAR.
3. Update `wdl/tasks/ProteinAnnotation/Foldseek.wdl` `default_docker`
   on the two ProstT5-touching tasks (or override via Terra method
   config) to `:10.0.2`. The Dockstore sync picks the new WDL up
   on the next push to `main`.
4. Resubmit the HvpViralProteinAnnotation workflow with call caching
   off (so the new image actually runs).
5. Read `stderr.log` from t_06 / t_07 — see "Reading the diagnostic
   output" below.

## CI semantics

This repo auto-bumps per-image VERSION on default-branch pushes
(`scripts/ci-bump-image-versions.sh`). **Do not edit the VERSION line
in `docker/foldseek-gpu/Makefile`** — any explicit edit is taken as
authoritative and bypasses the patch-bump (see `AGENTS.md` ► Auto-
versioning rules).

Procedure:

- Commit the patch + Dockerfile only. Leave VERSION at 10.0.1.
- Push the branch, open PR, merge.
- `docker.yml` detects the changed image, builds + scans, then
  `ci-bump-image-versions.sh` writes `VERSION = 10.0.2` and commits
  it back with `chore(docker): bump versions for foldseek-gpu [skip
  ci]`.
- The image is published as:
  - `ghcr.io/broadinstitute/hvp-lr/foldseek-gpu:10.0.2`
  - `us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.2`
- The `:latest` and `:main` tags also move.

## Wiring the WDL

`wdl/tasks/ProteinAnnotation/Foldseek.wdl` has two `default_docker`
inputs that need to advance:

- `FoldseekCreateDbFromFasta.default_docker` (t_06_FoldDb hits this)
- `FoldseekSearch.default_docker` (t_07_Search hits this through
  easy-search, which internally calls structcreatedb for the FASTA
  query)

Use the GAR tag — `ghcr.io` digests break Cromwell call caching on
GCPBATCH (see `MEMORY.md` ► HVP-LR ghcr.io callcache).

```diff
- default_docker = "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.1"
+ default_docker = "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.2"
```

If the WDL VERSION line in `Foldseek.wdl` (or any importing pipeline)
hasn't been touched in the PR yet, leave it — `cd.yml` handles the
repo-root VERSION + tag bump on merge.

### Option B: method config override (skip WDL edit)

For a one-shot validation run without re-syncing Dockstore:

```
HvpViralProteinAnnotation.FoldseekCreateDbFromFasta.docker_image
    = "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.2"
HvpViralProteinAnnotation.FoldseekSearch.docker_image
    = "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.2"
```

(Adjust input names if Foldseek.wdl exposes them under different
labels — the task surface defines `docker_image` via `default_docker`
through the standard HVP-LR runtime block).

## Resubmitting

Critical: **disable call caching** for the resubmit. Without it,
Cromwell happily replays the 10.0.1 cache hit and the new image
never runs. In Terra's submission dialog uncheck "Use call caching".

Use the padded BFVD target DB (already wired):
`gs://hvp-tech-dev/databases/foldseek/bfvd_foldseekdb_padded.tar.gz`.
The padded DB is independent of the verbose image — both fixes are
in flight together.

## Reading the diagnostic output

Stderr from t_06_FoldDb and t_07_Search will now include `[hvp-lr]`
prefixed lines. Grep for `[hvp-lr]` to isolate them. Expected
healthy output on a T4 host:

```
[hvp-lr] LlamaInitGuard: backend_init begin (verbose=1, quiet=0)
[hvp-lr] LlamaInitGuard: backend_init done
[hvp-lr] LlamaInitGuard: numa_init done (DISABLED)
[hvp-lr] LlamaInitGuard: ggml backend device count = 2
[hvp-lr] LlamaInitGuard:   dev[0] name=CUDA0 type=GPU desc=Tesla T4
[hvp-lr] LlamaInitGuard:   dev[1] name=CPU   type=CPU desc=...
[hvp-lr] ProstT5::getDevices: enumerating 2 backend device(s)
[hvp-lr] ProstT5::getDevices:   keep[0] name=CUDA0 type=GPU desc=Tesla T4
[hvp-lr] ProstT5::getDevices:   keep[1] name=CPU   type=CPU desc=...
[hvp-lr] ProstT5::getDevices: returning 2 device name(s)
[hvp-lr] ProstT5Model: ctor model_file=... device='CUDA0'
[hvp-lr] ProstT5Model: parse_device_list returned 2 entries (incl. nullptr terminator)
[hvp-lr] ProstT5Model:   parsed[0] name=CUDA0 type=GPU desc=Tesla T4
[hvp-lr] ProstT5Model:   parsed[1] = nullptr (terminator)
[hvp-lr] ProstT5Model: gpus=1 n_gpu_layers=24 use_mmap=1
[hvp-lr] ProstT5Model: llama_load_model_from_file begin
[llama:INFO] llm_load_tensors: offloaded 24/25 layers to GPU
[hvp-lr] ProstT5Model: llama_load_model_from_file OK
[hvp-lr] ProstT5: ctx ctor n_threads=N n_batch=2048 n_ubatch=2048 n_ctx=2048 embeddings=1 attn=NON_CAUSAL
[hvp-lr] ProstT5: llama_new_context_with_model OK
```

What the lines tell us:

| Symptom in stderr | Diagnosis |
|---|---|
| `dev[*] name=CPU` only — no CUDA device | GPU never attached. Check WDL `gpuType`, `gpuCount`, Cromwell GCPBATCH backend wiring. |
| CUDA device listed, but `parsed[0] = nullptr` only | `device` parameter empty or filtered out before parse — `--gpu 1` flag missing or stripped. |
| `gpus=1 n_gpu_layers=24` but no `[llama:INFO] llm_load_tensors: offloaded ... to GPU` lines | llama declined the GPU offload silently; very likely driver/CUDA mismatch (host driver too old for CUDA 12.6 runtime). |
| `llama_load_model_from_file FAILED` | Model file unreadable, corrupt, or wrong format. |
| Everything looks healthy but t_06 runtime still >> 10× T4 baseline | Throughput regression downstream of ctx creation — look at `n_ubatch / n_batch` and `[llama:*]` per-encode logs (each `predict()` call now triggers cached logging). |

## Verifying the build locally (optional)

```sh
cd docker/foldseek-gpu
make build_no_cache
docker run --rm --gpus all foldseek-gpu:10.0.1 structcreatedb \
    --help 2>&1 | head -20
```

(`structcreatedb --help` exits before ProstT5 init, so it won't
exercise the diagnostic lines. Run a real createdb against a tiny
FASTA + the ProstT5 weights to confirm `[hvp-lr]` lines appear.)

## Cleanup

The patch is a permanent diagnostic — leave it in. The `HVP_LR_QUIET=1`
env var suppresses llama internal logs if their volume becomes a
problem, but the `[hvp-lr]` markers will continue to print.
