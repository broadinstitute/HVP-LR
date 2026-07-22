# MMseqs2 built from source, pinned to upstream release tag 18-8cc5c
# (commit 8cc5ce367b5638c4306c2d7cfc652dd099a4643f).
#
# Multi-arch: linux/amd64 (SSE2 + SSE4.1 + AVX2 triple binaries + runtime
# dispatch wrapper) and linux/arm64 (single aarch64 binary, direct
# entrypoint). Runs natively on Linux/AMD64 and on Apple Silicon via Docker
# Desktop's linux/arm64 VM. No GPU/CUDA support — re-enable upstream's
# GPU=1 path in a sibling mmseqs2-gpu image if/when needed.
#
# Sibling to docker/foldseek/ — same Steinegger-lab build pattern. MMseqs2
# adds an SSE2 fallback variant (foldseek skips it). MMseqs2 has no Rust
# dependency, so the builder stage is simpler than foldseek's.
#
# Upstream Dockerfile reference:
#   https://github.com/soedinglab/MMseqs2/blob/<SHA>/Dockerfile
ARG MMSEQS_COMMIT=8cc5ce367b5638c4306c2d7cfc652dd099a4643f

ARG APP=mmseqs
FROM --platform=$BUILDPLATFORM debian:bookworm-slim AS source
ARG MMSEQS_COMMIT
RUN apt-get update \
    && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN git clone https://github.com/soedinglab/MMseqs2.git \
    && cd MMseqs2 \
    && git checkout "$MMSEQS_COMMIT" \
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

WORKDIR /opt/build
COPY --from=source /opt/MMseqs2/ .

# Three amd64 variants — sse2/sse41/avx2 — plus a single arm64 variant.
# Placeholder touches at the end of each branch keep the multi-binary COPY
# in the final stage architecture-agnostic.
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      mkdir -p build_$TARGETARCH/src; \
      cd /opt/build/build_$TARGETARCH; \
      CC=aarch64-linux-gnu-gcc CXX=aarch64-linux-gnu-g++ \
        cmake -DHAVE_ARM8=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 \
              -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j "$(nproc --all)"; \
      mv src/${APP} /opt/build/${APP}_arch; \
      touch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2; \
    else \
      mkdir -p build_sse2/src build_sse41/src build_avx2/src; \
      cd /opt/build/build_sse2; \
      cmake -DHAVE_SSE2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j "$(nproc --all)"; \
      mv src/${APP} /opt/build/${APP}_sse2; \
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
ARG MMSEQS_COMMIT

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/soedinglab/MMseqs2" \
      org.opencontainers.image.revision="${MMSEQS_COMMIT}" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.description="MMseqs2 sequence search and clustering, source-built and pinned" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /usr/local/bin/

# /usr/local/bin/mmseqs is both the user-callable command on PATH and the
# image's ENTRYPOINT. On amd64 it is the runtime-dispatch wrapper that picks
# mmseqs_sse2 vs mmseqs_sse41 vs mmseqs_avx2 based on /proc/cpuinfo; on
# arm64 it is a direct symlink to mmseqs_arch (the wrapper does not apply
# there).
COPY --from=source /opt/MMseqs2/util/${APP}_wrapper.sh /usr/local/bin/${APP}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      gawk bash grep procps libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar aria2 ca-certificates \
    && rm -rf /var/lib/apt/lists/*; \
    chmod +x /usr/local/bin/${APP}; \
    if [ "$TARGETARCH" = "arm64" ]; then \
      rm -f /usr/local/bin/${APP}; \
      ln -s /usr/local/bin/${APP}_arch /usr/local/bin/${APP}; \
    fi

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image. MMseqs2
# itself is GPL-3.0-or-later — see image label above.
COPY LICENSE /opt/hvp-lr/LICENSE

ENTRYPOINT ["/usr/local/bin/mmseqs"]
