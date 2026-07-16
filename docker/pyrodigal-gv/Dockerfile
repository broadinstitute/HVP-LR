# Pyrodigal-gv — Python/Cython wrapper around Prodigal-gv (Prodigal fork
# tuned for giant viruses + jumbo phages, includes alt-code training models
# for NCLDV, Topaz/Agate code-15, gut phages). Used for ORF calling on viral
# nucleotide sequences where standard Prodigal would miss 10-30% of genes
# (TAG-as-Q, alt start codons).
#
# Multi-arch: linux/amd64 + linux/arm64. Pyrodigal ships manylinux2_28 wheels
# for both architectures so no Cython rebuild needed in this image.
#
# Pinned to pyrodigal-gv 0.3.2 (released 2024-08-15) which depends on
# pyrodigal ~=3.0. We pin pyrodigal exactly to 3.7.1 to keep the runtime
# deterministic. Bump VERSION in the Makefile and both pins below together.
ARG PYRODIGAL_GV_VERSION=0.3.2
ARG PYRODIGAL_VERSION=3.7.1

FROM --platform=$TARGETPLATFORM python:3.12-slim-bookworm AS runtime
ARG PYRODIGAL_GV_VERSION
ARG PYRODIGAL_VERSION

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/althonos/pyrodigal-gv" \
      org.opencontainers.image.version="${PYRODIGAL_GV_VERSION}" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.description="Pyrodigal-gv viral ORF caller (Prodigal-gv via Cython)" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

# Runtime utilities (procps for `free -g` in WDL preamble; standard data-handling tools).
# Wheels are pulled directly; no compiler toolchain needed in the runtime image.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      bash grep gawk procps wget tar gzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# pyrodigal ~3.7 wheels carry triple-tagged manylinux2014 + 2_17 + 2_28
# compatibility for both x86_64 and aarch64. python:3.12-slim is glibc-based,
# so pip resolves the manylinux wheel on both architectures without building
# from source.
RUN pip install --no-cache-dir \
      "pyrodigal==${PYRODIGAL_VERSION}" \
      "pyrodigal-gv==${PYRODIGAL_GV_VERSION}"

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image. Pyrodigal-gv
# itself is GPL-3.0-or-later — see image label above.
COPY LICENSE /opt/hvp-lr/LICENSE

ENTRYPOINT ["pyrodigal-gv"]
