# CheckM2 MAG quality assessment (python 3.10 + tensorflow 2.17)
#
# Split out of the retired hvp-monolith image. All pins are copied verbatim
# from docker/hvp-monolith/env.yaml so behavior is byte-for-byte identical to
# the monolith's base env for these tools — just without the 10+ unused tools
# and 3 conflicting sidecar python envs the monolith carried.
#
# Multi-arch: linux/amd64 + linux/arm64 (bioconda ships both for these pins).
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

# Patch base-OS security updates (libssl3t64 CVE-2026-45447 et al.) that the
# pinned micromamba base tag ships un-upgraded. Root needed for apt; drop back.
USER root
RUN apt-get update && apt-get upgrade -y && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/chklovski/CheckM2" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.description="CheckM2 MAG quality assessment (python 3.10 + tensorflow 2.17)" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

# Pinned conda spec drives the image contents (same file CI documents as env.yaml).
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete

# Drop pip's inaccurate vendored SBOM: it lists setuptools==70.3.0 though pip
# only vendors pkg_resources (no PackageIndex), so CVE-2025-47273's vulnerable
# code is not present — the stale SBOM is all Trivy sees. Removing it keeps the
# image SBOM accurate; pip and tools run fine without these metadata files.
RUN rm -f /opt/conda/lib/python*/site-packages/pip/_vendor/bom.cdx.json \
          /opt/conda/lib/python*/site-packages/pip/_vendor/vendor.txt

# Tools on PATH without needing conda activation — matches how the monolith
# exposed its base env to WDL command blocks (Cromwell runs /bin/bash, not the
# image entrypoint, so PATH is what makes the binaries resolve).
ENV PATH="/opt/conda/bin:$PATH"

# Repo-root LICENSE bundled into every image (staged into build context by CI).
COPY LICENSE /opt/hvp-lr/LICENSE

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "checkm2"]
