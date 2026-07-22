// HVP-LR batched-encode verifier for ProstT5 / llama.cpp t5encoder.
//
// Pre-flight check for the (proposed) batched-encode patch: validates
// that llama.cpp's t5encoder NON_CAUSAL attention handles multi-seq
// batches correctly on the active backend (target: CUDA). If batched
// per-token embeddings match single-seq per-token embeddings within
// float-noise tolerance for every seq in a K-packed batch, then the
// seq_id-based attention mask is masking cross-seq tokens correctly and
// it is safe to pursue the batched-encode patch. If they diverge, the
// batched-encode path would silently produce wrong 3Di strings.
//
// Usage:
//   verify_batched_encode <model.gguf> [K] [L]
//     K  — number of sequences to pack (default 4)
//     L  — token length per sequence (default 64)
//
// Exit code: 0 on PASS, 1 on FAIL, 2 on argument / model load error.

#include "llama.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static void batch_add(llama_batch & batch,
                      llama_token   tok,
                      llama_pos     pos,
                      llama_seq_id  seq_id) {
    int32_t idx = batch.n_tokens;
    batch.token[idx]    = tok;
    batch.pos[idx]      = pos;
    batch.n_seq_id[idx] = 1;
    batch.seq_id[idx][0] = seq_id;
    batch.logits[idx]   = 1;
    batch.n_tokens++;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <model.gguf> [K] [L]\n", argv[0]);
        return 2;
    }
    const char * modelPath = argv[1];
    int K = (argc > 2) ? std::atoi(argv[2]) : 4;
    int L = (argc > 3) ? std::atoi(argv[3]) : 64;
    if (K < 2 || K > 32 || L < 4 || L > 512) {
        std::fprintf(stderr, "K must be in [2,32], L in [4,512] (got K=%d L=%d)\n", K, L);
        return 2;
    }

    llama_backend_init();

    auto mparams = llama_model_default_params();
    mparams.n_gpu_layers = 999;
    llama_model * model = llama_load_model_from_file(modelPath, mparams);
    if (model == nullptr) {
        std::fprintf(stderr, "ERROR: failed to load model from %s\n", modelPath);
        llama_backend_free();
        return 2;
    }

    auto cparams = llama_context_default_params();
    cparams.n_ctx          = 2048;
    cparams.n_batch        = 2048;
    cparams.n_ubatch       = 2048;
    cparams.embeddings     = true;
    cparams.attention_type = LLAMA_ATTENTION_TYPE_NON_CAUSAL;
    // Match the runtime ProstT5 context config: foldseek's ProstT5.cpp sets
    // flash_attn=true on every encode context (multi-context-prostt5.patch).
    // FA kernels and dense-matmul kernels can mask differently — the
    // verifier MUST exercise the same path foldseek will use in production.
    // Override at runtime with HVP_LR_VERIFY_FLASH_ATTN=0 to test both paths.
    const char * faenv = std::getenv("HVP_LR_VERIFY_FLASH_ATTN");
    cparams.flash_attn = (faenv == nullptr) ? true : (std::atoi(faenv) != 0);
    std::fprintf(stderr, "[verify] flash_attn=%d\n", cparams.flash_attn ? 1 : 0);

    llama_context * ctx = llama_new_context_with_model(model, cparams);
    if (ctx == nullptr) {
        std::fprintf(stderr, "ERROR: failed to create context\n");
        llama_free_model(model);
        llama_backend_free();
        return 2;
    }

    int n_embd = llama_n_embd(model);
    std::fprintf(stderr, "[verify] model=%s K=%d L=%d n_embd=%d\n",
                 modelPath, K, L, n_embd);

    // Build K distinct token sequences. Stay clear of low-numbered control
    // tokens (0-3 are typically reserved). Token id 4..127 should be safe
    // for the ProstT5 SentencePiece vocab (AA tokens + 3Di tokens).
    std::vector<std::vector<llama_token>> seqs(K);
    for (int k = 0; k < K; ++k) {
        seqs[k].reserve(L);
        for (int i = 0; i < L; ++i) {
            llama_token t = 4 + ((k * 7 + i * 3) % 64);
            seqs[k].push_back(t);
        }
    }

    // --- Path A: single-seq encode per sequence ---
    // No llama_kv_cache_clear between calls: the pinned llama.cpp guards
    // the public llama_kv_cache_* functions on llama_context behind #if 0
    // (lines 21487-22609 of lib/prostt5/src/llama.cpp). Safe to skip — the
    // t5encoder NON_CAUSAL graph is per-pass scratch, not a growing cache;
    // foldseek's ProstT5::predict (src/strucclustutils/ProstT5.cpp:64)
    // already loops llama_encode per sequence with no clear between.
    std::vector<std::vector<float>> singleEmbd(K);
    for (int k = 0; k < K; ++k) {
        llama_batch b = llama_batch_init(L, 0, 1);
        for (int i = 0; i < L; ++i) {
            batch_add(b, seqs[k][i], (llama_pos) i, (llama_seq_id) 0);
        }
        if (llama_encode(ctx, b) < 0) {
            std::fprintf(stderr, "ERROR: llama_encode failed in single-seq path (k=%d)\n", k);
            llama_batch_free(b);
            llama_free(ctx);
            llama_free_model(model);
            llama_backend_free();
            return 1;
        }
        singleEmbd[k].resize((size_t) L * n_embd);
        for (int i = 0; i < L; ++i) {
            float * e = llama_get_embeddings_ith(ctx, i);
            if (e == nullptr) {
                std::fprintf(stderr, "ERROR: get_embeddings_ith returned null (single k=%d i=%d)\n", k, i);
                llama_batch_free(b);
                llama_free(ctx);
                llama_free_model(model);
                llama_backend_free();
                return 1;
            }
            std::memcpy(&singleEmbd[k][(size_t) i * n_embd], e, n_embd * sizeof(float));
        }
        llama_batch_free(b);
    }

    // --- Path B: packed K-seq encode in one batch ---
    llama_batch packed = llama_batch_init(K * L, 0, K);
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < L; ++i) {
            batch_add(packed, seqs[k][i], (llama_pos) i, (llama_seq_id) k);
        }
    }
    if (llama_encode(ctx, packed) < 0) {
        std::fprintf(stderr, "ERROR: llama_encode failed in packed path (K=%d L=%d)\n", K, L);
        llama_batch_free(packed);
        llama_free(ctx);
        llama_free_model(model);
        llama_backend_free();
        return 1;
    }

    std::vector<std::vector<float>> packedEmbd(K);
    for (int k = 0; k < K; ++k) {
        packedEmbd[k].resize((size_t) L * n_embd);
    }
    for (int p = 0; p < K * L; ++p) {
        int k = p / L;
        int i = p % L;
        float * e = llama_get_embeddings_ith(ctx, p);
        if (e == nullptr) {
            std::fprintf(stderr, "ERROR: get_embeddings_ith returned null (packed p=%d)\n", p);
            llama_batch_free(packed);
            llama_free(ctx);
            llama_free_model(model);
            llama_backend_free();
            return 1;
        }
        std::memcpy(&packedEmbd[k][(size_t) i * n_embd], e, n_embd * sizeof(float));
    }
    llama_batch_free(packed);

    // --- Compare ---
    double overallMax  = 0.0;
    double overallSum  = 0.0;
    long   overallN    = 0;
    for (int k = 0; k < K; ++k) {
        double seqMax = 0.0;
        double seqSum = 0.0;
        long   seqN   = 0;
        for (size_t j = 0; j < (size_t) L * n_embd; ++j) {
            double d = std::fabs((double) singleEmbd[k][j] - (double) packedEmbd[k][j]);
            if (d > seqMax) seqMax = d;
            seqSum += d;
            ++seqN;
        }
        std::fprintf(stderr, "[verify] seq[%d] max_abs_diff=%.6g mean_abs_diff=%.6g\n",
                     k, seqMax, seqSum / (double) seqN);
        if (seqMax > overallMax) overallMax = seqMax;
        overallSum += seqSum;
        overallN   += seqN;
    }
    std::fprintf(stderr, "[verify] overall max_abs_diff=%.6g mean_abs_diff=%.6g\n",
                 overallMax, overallSum / (double) overallN);

    // Tolerance: float32 noise from non-deterministic CUDA matmul reduction
    // order is well under 1e-3 for embeddings of this scale. A correctly-
    // masked batched path should be near-bit-identical; a broken mask
    // produces large divergence (cross-seq attention leaks into logits).
    const double tol = 1e-3;
    bool pass = overallMax < tol;
    std::fprintf(stderr, "[verify] %s — max_abs_diff %s tol=%.0e\n",
                 pass ? "PASS" : "FAIL",
                 pass ? "<" : ">=", tol);

    llama_free(ctx);
    llama_free_model(model);
    llama_backend_free();
    return pass ? 0 : 1;
}
