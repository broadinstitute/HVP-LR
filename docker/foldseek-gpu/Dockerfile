# Foldseek built from source with CUDA enabled, pinned to upstream master
# @ 718d4217. Sibling to docker/foldseek/ but adds NVIDIA CUDA 12.6 build
# deps so foldseek's ProstT5 inference path can use GPUs when present.
#
# Architecture: linux/amd64 ONLY. CUDA has no arm64 build for the bookworm
# repo we target.
#
# Compute capabilities baked in (foldseek upstream selection):
#   75-real;80-real;86-real;89-real;90
# i.e. Turing (T4), Ampere (A100, RTX 30xx), Ada Lovelace (RTX 40xx, L4),
# Hopper (H100). Older cards (Pascal P100, Volta V100) will not load this
# binary — use the CPU image (docker/foldseek/) for those.
#
# Runtime requirement: host must have NVIDIA drivers + nvidia-container-toolkit
# and the container must be launched with --gpus all (or equivalent).
# When run on a host without GPU access, foldseek's GPU code paths will fail
# at import; CPU-only callers should use the sibling docker/foldseek/ image.
ARG FOLDSEEK_COMMIT=718d42176d2f67d36a60866fedfb881f8d5a7ebf

ARG APP=foldseek
FROM --platform=$BUILDPLATFORM debian:bookworm-slim AS source
ARG FOLDSEEK_COMMIT
RUN apt-get update \
    && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN git clone https://github.com/steineggerlab/foldseek.git \
    && cd foldseek \
    && git checkout "$FOLDSEEK_COMMIT" \
    && git submodule update --init --recursive

FROM --platform=$BUILDPLATFORM debian:bookworm-slim AS builder
ARG TARGETARCH
ARG APP

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      build-essential curl xxd cmake ninja-build ca-certificates patch \
      zlib1g-dev libbz2-dev libatomic1 \
    && rm -rf /var/lib/apt/lists/*

# NVIDIA CUDA toolkit (matches upstream's GPU=1 build path).
RUN curl -fsSL -o /tmp/cuda-keyring.deb \
      https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/cuda-keyring_1.1-1_all.deb \
    && dpkg -i /tmp/cuda-keyring.deb \
    && rm /tmp/cuda-keyring.deb \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      cuda-nvcc-12-6 cuda-cudart-dev-12-6 libcublas-dev-12-6 \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/usr/local/cuda-12.6/bin:${PATH}"

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
      sh -s -- --profile minimal --default-host x86_64-unknown-linux-gnu -y
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /opt/build
COPY --from=source /opt/foldseek/ .

# HVP-LR diagnostic patch: forces verbose llama init and adds [hvp-lr]
# traces around ProstT5 device/model/context construction. Required to
# distinguish "GPU attached but inference fell back to CPU" from "GPU
# never attached" in production stderr — upstream silences every llama
# init signal because the MMseqs2 verbosity cap (3) is one below the
# threshold that toggles llama logs.
COPY patches/ /opt/patches/
RUN patch -p1 < /opt/patches/verbose-prostt5.patch

# HVP-LR multi-context GPU concurrency patch. Hoists ProstT5Model out of
# the OMP parallel region, runs N llama_contexts in parallel on a shared
# (mmap'd) model. Defaults to N=8; override at runtime via env
# HVP_LR_PROSTT5_GPU_CTXS=<int>. See patches/README.md for VRAM math.
RUN patch -p1 < /opt/patches/multi-context-prostt5.patch

# HVP-LR batched-encode pre-flight verifier. Standalone C++ binary that
# validates llama.cpp t5encoder NON_CAUSAL multi-seq attention masking on
# the active backend (target: CUDA). Required gate before pursuing the
# batched-encode patch — if the mask is broken, batched encode would
# silently produce wrong 3Di strings. See verify/README.md.
COPY verify/ /opt/verify/
RUN cp /opt/verify/verify_batched_encode.cpp src/strucclustutils/ \
    && patch -p1 < /opt/patches/verify-batched-encode.patch

# Single AVX2 + CUDA build (upstream's GPU=1 path). The wrapper script is not
# used; the entrypoint is symlinked directly to foldseek_avx2 in the final
# stage. SSE4.1 + arch placeholders are touched so the multi-binary COPY in
# the final stage succeeds.
# && (not ;) so any failure in cmake/ninja/mv aborts the layer instead of
# silently producing a partial image that breaks the final-stage COPY.
RUN mkdir -p build_cuda/src \
    && cd /opt/build/build_cuda \
    && cmake -GNinja \
          -DHAVE_AVX2=1 \
          -DENABLE_CUDA=1 \
          -DCMAKE_CUDA_ARCHITECTURES="75-real;80-real;86-real;89-real;90" \
          -DHAVE_MPI=0 -DHAVE_TESTS=0 \
          -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. .. \
    && ninja -j "$(nproc --all)" \
    && mv src/${APP} /opt/build/${APP}_avx2 \
    && mv src/verify_batched_encode /opt/build/verify_batched_encode \
    && touch /opt/build/${APP}_arch /opt/build/${APP}_sse41

FROM debian:bookworm-slim
ARG APP
ARG FOLDSEEK_COMMIT

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/steineggerlab/foldseek" \
      org.opencontainers.image.revision="${FOLDSEEK_COMMIT}" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.description="Foldseek structural aligner with CUDA-enabled ProstT5 inference, source-built and pinned" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

# CUDA runtime libs (libcudart, libcublas) so the binary loads on GPU hosts
# without needing the host to mount its own CUDA runtime in.
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl ca-certificates \
    && curl -fsSL -o /tmp/cuda-keyring.deb \
      https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/cuda-keyring_1.1-1_all.deb \
    && dpkg -i /tmp/cuda-keyring.deb \
    && rm /tmp/cuda-keyring.deb \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      cuda-cudart-12-6 libcublas-12-6 \
      gawk bash grep procps libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar aria2 \
    && apt-get purge -y curl \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
ENV LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /opt/build/verify_batched_encode /usr/local/bin/

# /usr/local/bin/foldseek is the user-callable command on PATH and the image
# ENTRYPOINT. CUDA build emits only the AVX2 variant; direct symlink (no
# wrapper-based runtime dispatch needed because CUDA hosts are always AVX2).
RUN ln -s /usr/local/bin/${APP}_avx2 /usr/local/bin/${APP}

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image. Foldseek
# itself is GPL-3.0-or-later — see image label above.
COPY LICENSE /opt/hvp-lr/LICENSE

ENTRYPOINT ["/usr/local/bin/foldseek"]
