# Foldseek built from source, pinned to upstream master @ 718d4217.
#
# Multi-arch: linux/amd64 (SSE4.1 + AVX2 dual binaries + runtime dispatch
# wrapper) and linux/arm64 (single aarch64 binary, direct entrypoint).
# Runs natively on Linux/AMD64 and on Apple Silicon via Docker Desktop's
# linux/arm64 VM. No GPU/CUDA support in this image — strip that path from
# upstream's Dockerfile; re-add if/when a GPU runner is added to CI.
#
# Upstream Dockerfile reference:
#   https://github.com/steineggerlab/foldseek/blob/<SHA>/Dockerfile
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

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      build-essential curl xxd cmake ca-certificates \
      zlib1g-dev libbz2-dev libatomic1 \
      crossbuild-essential-$TARGETARCH zlib1g-dev:$TARGETARCH libbz2-dev:$TARGETARCH \
    && rm -rf /var/lib/apt/lists/*

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      ARCH=aarch64; \
    else \
      ARCH=x86_64; \
    fi; \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
      sh -s -- --profile minimal --default-host ${ARCH}-unknown-linux-gnu -y
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /opt/build
COPY --from=source /opt/foldseek/ .

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      mkdir -p build_$TARGETARCH/src; \
      cd /opt/build/build_$TARGETARCH; \
      CC=aarch64-linux-gnu-gcc CXX=aarch64-linux-gnu-g++ \
        cmake -DHAVE_ARM8=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 \
              -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j "$(nproc --all)"; \
      mv src/${APP} /opt/build/${APP}_arch; \
      touch /opt/build/${APP}_sse41 /opt/build/${APP}_avx2; \
    else \
      mkdir -p build_sse41/src build_avx2/src; \
      cd /opt/build/build_sse41; \
      cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j "$(nproc --all)"; \
      mv src/${APP} /opt/build/${APP}_sse41; \
      cd /opt/build/build_avx2; \
      cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j "$(nproc --all)"; \
      mv src/${APP} /opt/build/${APP}_avx2; \
      touch /opt/build/${APP}_arch; \
    fi

FROM debian:bookworm-slim
ARG TARGETARCH
ARG APP
ARG FOLDSEEK_COMMIT

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/steineggerlab/foldseek" \
      org.opencontainers.image.revision="${FOLDSEEK_COMMIT}" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.description="Foldseek structural aligner, source-built and pinned" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /usr/local/bin/

# /usr/local/bin/foldseek is both the user-callable command on PATH and the
# image's ENTRYPOINT. On amd64 it is the runtime-dispatch wrapper that picks
# foldseek_sse41 vs foldseek_avx2 based on /proc/cpuinfo; on arm64 it is a
# direct symlink to foldseek_arch (the wrapper does not apply there).
COPY --from=source /opt/foldseek/util/${APP}_wrapper.sh /usr/local/bin/${APP}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      gawk bash grep procps libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar aria2 ca-certificates \
    && rm -rf /var/lib/apt/lists/*; \
    if [ "$TARGETARCH" = "arm64" ]; then \
      rm -f /usr/local/bin/${APP}; \
      ln -s /usr/local/bin/${APP}_arch /usr/local/bin/${APP}; \
    fi

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image. Foldseek
# itself is GPL-3.0-or-later — see image label above.
COPY LICENSE /opt/hvp-lr/LICENSE

ENTRYPOINT ["/usr/local/bin/foldseek"]
