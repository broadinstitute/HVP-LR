# HVP-LR batched-encode verifier

Pre-flight gate for the (proposed) ProstT5 batched-encode patch. Validates
that llama.cpp's t5encoder NON_CAUSAL attention handles a multi-seq batch
correctly on the active backend (target: CUDA).

The verifier compares two paths against the same model + same tokens:

- **Path A** — K single-seq `llama_encode` calls (the current foldseek
  inference path), KV cache cleared between.
- **Path B** — one packed K-seq `llama_encode` call with `seq_id={0..K-1}`
  in a single `llama_batch`.

If per-token embeddings match within float-noise tolerance (`max_abs_diff
< 1e-3`) for every seq, the seq_id-based attention mask is doing its job
and the batched-encode patch is safe to pursue. If they diverge, the
batched-encode path would silently produce wrong 3Di strings.

## Building

The verifier is built as a CMake target inside the foldseek-gpu image
(see `patches/verify-batched-encode.patch`). The binary lands at
`/usr/local/bin/verify_batched_encode` in the final image.

## Running

Requires a GPU host (CUDA backend behavior is what's being validated; the
CPU backend may mask correctly even if CUDA doesn't). Run against the
ProstT5 GGUF that production uses:

```sh
docker run --rm --gpus all \
    -v "$PWD":/work -w /work \
    --entrypoint /usr/local/bin/verify_batched_encode \
    foldseek-gpu:10.0.2 \
    prostt5-f16.gguf
```

Defaults: `K=4` sequences of `L=64` tokens. Override:

```sh
docker run ... verify_batched_encode prostt5-f16.gguf 8 128
```

Bounds: `K in [2,32]`, `L in [4,512]`. The packed path consumes ~K*L
KV-cache slots and one compute buffer sized for K*L; the container
defaults (`n_ctx=n_batch=n_ubatch=2048`) accommodate up to K*L=2048.

## Exit codes

| Code | Meaning |
|------|---------|
| 0    | PASS — batched encode is bit-equivalent (within tol) to single-seq |
| 1    | FAIL — divergence > tol, batched-encode patch NOT safe |
| 2    | Argument or model-load error |

## Output

Per-seq and overall `max_abs_diff` / `mean_abs_diff` lines on stderr,
prefixed `[verify]`. The PASS/FAIL verdict is the last `[verify]` line.

## Interpreting FAIL

If `max_abs_diff` is in the `1e-2` to `1e+0` range, cross-seq attention
is leaking: the NON_CAUSAL mask in `llama.cpp/src/llama-graph.cpp`
(specifically the encoder branch) is not partitioning tokens by
`seq_id`. Do NOT pursue the batched-encode patch on this CUDA backend
without first patching the mask in the vendored llama.cpp tree.
