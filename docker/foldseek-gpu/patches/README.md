# foldseek-gpu HVP-LR diagnostic patches

`patch -p1` files applied to the foldseek source tree during the
`docker/foldseek-gpu/` builder stage (see Dockerfile around line 60).
The patches run after `COPY --from=source /opt/foldseek/ .` and before
the cmake invocation, so the modified `.cpp` files compile into the
shipped `foldseek_avx2` binary.

Patches are tied to the upstream commit pinned in `Dockerfile`'s
`FOLDSEEK_COMMIT` ARG (currently
`718d42176d2f67d36a60866fedfb881f8d5a7ebf`). Bumping the foldseek pin
without regenerating the patches will fail the `RUN patch -p1` step.

## verbose-prostt5.patch

Forces verbose llama / ggml init inside the ProstT5 inference path and
adds `[hvp-lr]` diagnostic traces around device enumeration, model
load, and context creation. Required to distinguish "GPU attached but
inference fell back to CPU" from "GPU never attached" in production
stderr.

**Why this is needed.** The upstream `LlamaInitGuard` (in
`src/strucclustutils/ProstT5.cpp`) installs a no-op llama log callback
unless its `verbose` parameter is true. The only caller in
`src/strucclustutils/structcreatedb.cpp` sets `verbose = (par.verbosity
> 3)`, but the MMseqs2 PARAM_V regex caps verbosity at 3
(`^[0-3]{1}$`), so the gate is unreachable from the CLI. Every llama
init signal — CUDA backend init failures, n_gpu_layers actually
applied, ggml device enumeration — is silently dropped.

**What the patch changes.**

1. `LlamaInitGuard` ctor passes llama logs through to stderr unless
   `HVP_LR_QUIET=1` is set (override available for the cases where
   we want clean logs). Emits an `[hvp-lr]` enumeration of every
   ggml backend device (name, type tag, description) immediately
   after `llama_backend_init()`.
2. `ProstT5Model` ctor logs the model file path, the resolved device
   string, every entry returned by `parse_device_list()`, the GPU
   count, the `n_gpu_layers` value finally written to mparams, and
   whether `llama_load_model_from_file()` returned a non-null model.
3. `ProstT5` (context) ctor logs `n_threads / n_batch / n_ubatch /
   n_ctx` and whether `llama_new_context_with_model()` returned a
   non-null context.
4. `ProstT5::getDevices()` logs every device it enumerates, including
   the type tag (GPU / CPU / ACCEL) and whether it was kept or
   skipped.
5. `structcreatedb.cpp` line 776 forces
   `LlamaInitGuard guard(true)` (the only place the upstream cap on
   `par.verbosity > 3` was checked).

All diagnostic lines are prefixed `[hvp-lr]` so they grep cleanly out
of mixed stderr.

## multi-context-prostt5.patch

GPU performance bundle for the ProstT5 inference path. Three orthogonal
changes that compose:

1. **Multi-context concurrency** — N `llama_context`s in parallel on one
   mmap'd `ProstT5Model` (default 8, env `HVP_LR_PROSTT5_GPU_CTXS`).
2. **Flash attention** — `cparams.flash_attn = true` so backends that
   ship a flash-attn kernel for the t5encoder NON_CAUSAL graph use it
   (CPU/CUDA fall back silently when unsupported).
3. **Dynamic n_ubatch** — picks the largest safe `cparams.n_ubatch`
   (default 2048, env `HVP_LR_PROSTT5_GPU_UBATCH`) from the free VRAM
   budget at startup. Smaller n_ubatch shrinks the per-ctx CUDA0
   compute buffer ~linearly without touching max-seq capacity
   (n_ctx + n_batch stay at 2048 so the encoder still accepts up-to-2046
   AA chunks paired with the upstream `prostt5SplitLength` default).

All three are bit-equivalent at the 3Di-string level — every sequence
still hits the same single-seq `llama_encode` path.

**Why this is needed.** Upstream's GPU path sets
`localThreads = devices.size()` (almost always 1 CUDA device on Terra
hosts) and constructs one `ProstT5Model` + `ProstT5` ctx per OMP
thread. Result: the whole ProstT5 inference loop runs single-threaded
sequentially at ~175 ms/seq on Tesla T4 even though the kernels are
GPU-resident. A 6.5k-ORF DB takes ~20 min where the underlying CUDA
throughput would support ~3 min at 8x concurrency.

**What the patch changes** (`src/strucclustutils/structcreatedb.cpp`,
GPU branch only):

1. Hoists `ProstT5Model` out of the `#pragma omp parallel` region.
   One mmap'd copy of the weights, one GPU upload, shared by every
   context.
2. Replaces `localThreads = devices.size()` with
   `localThreads = min(HVP_LR_PROSTT5_GPU_CTXS, par.threads)`,
   default 8.
3. Each OMP thread constructs its own `ProstT5(model, 1)` ctx (own
   KV cache + compute buffer), dispatched across the ORF DB via the
   existing `#pragma omp for schedule(dynamic, 1)`.
4. CPU path (`par.gpu == 0`) untouched — still uses the per-thread
   `ProstT5Model` construction path, since the fork-per-process runner
   in `ProstT5ForkRunner.h` already handles CPU parallelism there.

**VRAM budget (measured, not estimated).** Per-context cost at
`n_ctx = n_batch = n_ubatch = 2048` (upstream defaults), F16 KV cache,
24 layers offloaded, observed in createdb.log on Quadro RTX 8000:

| Component        | Per-ctx | Shared |
|------------------|---------|--------|
| Model weights (`CUDA0 model buffer size`) |         | ~2.3 GB |
| CUDA runtime + cuBLAS + foldseek baseline |         | ~0.6 GB |
| KV cache (`n_embd_k_gqa=4096`, 24 layers, F16) | ~1.5 GB | |
| Compute buffer (`CUDA0 compute buffer size` at n_ubatch=2048) | ~2.26 GB | |
| **Per-ctx total**                        | **~3.8 GB** | |

That's ~15x what the initial design estimate assumed. Per-GPU ceiling
(VRAM-bound, before SM-occupancy saturation):

| GPU | VRAM | Max N @ n_ubatch=2048 | Recommended N | n_ubatch fallback |
|-----|------|-----------------------|---------------|-------------------|
| Tesla T4 | 16 GB | 3 | 3 | shrink ubatch → ~1024 → N=5-6 |
| L4 | 24 GB | 5 | 5 | shrink ubatch → ~1536 → N=8 |
| A10G | 24 GB | 5 | 5 | as L4 |
| Quadro RTX 8000 | 48 GB | 12 | 8 | full 2048 |
| A100-40 | 40 GB | 9 | 8 | full 2048 |
| A100-80 | 80 GB | 20 | 12 | full 2048 |

**Overrides.**
- `HVP_LR_PROSTT5_GPU_CTXS=<int>` — number of contexts (default 8).
- `HVP_LR_PROSTT5_GPU_UBATCH=<int>` in `[256, 2048]` — pins physical
  sub-batch size. If unset, the patch queries free VRAM at startup
  (via `ggml_backend_dev_memory()` over the selected device) and
  snaps to the largest safe step in `{2048, 1536, 1024, 768, 512, 256}`
  given the chosen N and the per-ctx KV + shared overhead budget.
  See `[hvp-lr] ProstT5 n_ubatch (auto|env|fallback)` log line.

**Expected speedup.** Sub-linear with N (memory-bandwidth contention
on the encoder weights). Target T4 N=3 + flash-attn → ~3x baseline.
L4 N=8 + flash-attn → ~5-6x. Numbers require empirical sweep on each
GPU class (the 36-seq dev test in the original validation was too
small to amortize per-ctx startup; need ≥10k seqs for meaningful
timing).

## Regenerating the patch

Patches are tracked as text. To edit:

```sh
cd /tmp
git clone https://github.com/steineggerlab/foldseek.git fs_src
cd fs_src
git checkout 718d42176d2f67d36a60866fedfb881f8d5a7ebf

# Snapshot then edit the files in place
cp src/strucclustutils/ProstT5.cpp        ../ProstT5.orig.cpp
cp src/strucclustutils/structcreatedb.cpp ../structcreatedb.orig.cpp
$EDITOR src/strucclustutils/ProstT5.cpp src/strucclustutils/structcreatedb.cpp

# Regenerate combined patch
diff -u --label a/src/strucclustutils/ProstT5.cpp \
        --label b/src/strucclustutils/ProstT5.cpp \
        ../ProstT5.orig.cpp src/strucclustutils/ProstT5.cpp \
  > <repo>/docker/foldseek-gpu/patches/verbose-prostt5.patch
diff -u --label a/src/strucclustutils/structcreatedb.cpp \
        --label b/src/strucclustutils/structcreatedb.cpp \
        ../structcreatedb.orig.cpp src/strucclustutils/structcreatedb.cpp \
  >> <repo>/docker/foldseek-gpu/patches/verbose-prostt5.patch

# Validate
patch --dry-run -p1 -d /tmp/fs_src \
    < <repo>/docker/foldseek-gpu/patches/verbose-prostt5.patch
```
