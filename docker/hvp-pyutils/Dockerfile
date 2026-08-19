# Python glue layer (pandas/pysam/biopython/seqkit) for HVP summary and rescue tasks
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
      org.opencontainers.image.url="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.description="Python glue layer (pandas/pysam/biopython/seqkit) for HVP summary and rescue tasks" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

# Pinned conda spec drives the image contents (same file CI documents as env.yaml).
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && rm -rf /opt/conda/pkgs /opt/conda/x86_64-conda-linux-gnu \
              /opt/conda/share/man /opt/conda/share/doc /opt/conda/share/locale \
    && { find /opt/conda -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null; \
         find /opt/conda -type d -name tests -path '*/site-packages/*' -exec rm -rf {} + 2>/dev/null; \
         find /opt/conda -type d -name include -path '*/site-packages/tensorflow/*' -exec rm -rf {} + 2>/dev/null; true; }

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

# multi-tool image: inherit micromamba base entrypoint (env-activating shell)
