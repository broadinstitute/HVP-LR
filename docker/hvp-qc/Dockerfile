# fastqc + multiqc — read/report QC tools carried by the retired hvp-monolith
# but not used by the active WDL pipelines. Kept available in a small own image.
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

# Base-OS security patch (libssl3t64 CVE-2026-45447 et al.); root for apt.
USER root
RUN apt-get update && apt-get upgrade -y && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/MultiQC/MultiQC" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.description="FastQC + MultiQC sequence QC reporting" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

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

# Drop pip's inaccurate vendored SBOM (see CVE-2025-47273 note in sibling images).
RUN rm -f /opt/conda/lib/python*/site-packages/pip/_vendor/bom.cdx.json \
          /opt/conda/lib/python*/site-packages/pip/_vendor/vendor.txt

ENV PATH="/opt/conda/bin:$PATH"
COPY LICENSE /opt/hvp-lr/LICENSE
# multi-tool image: inherit micromamba base entrypoint
