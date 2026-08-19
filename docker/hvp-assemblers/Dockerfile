# flye + hifiasm_meta + metamdbg — alternate metagenome assemblers carried by
# the monolith but not used by the active myloasm-based pipeline.
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

USER root
RUN apt-get update && apt-get upgrade -y && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/fenderglass/Flye" \
      org.opencontainers.image.licenses="BSD-3-Clause / MIT (per tool)" \
      org.opencontainers.image.description="Flye + hifiasm_meta + metaMDBG metagenome assemblers" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete

# metamdbg per-arch: bioconda ships 1.4 on linux-64, 1.2 on linux-aarch64.
ARG TARGETARCH
RUN set -e; \
    case "$TARGETARCH" in \
        amd64) METAMDBG_VER=1.4 ;; \
        arm64) METAMDBG_VER=1.2 ;; \
        *)     METAMDBG_VER=1.4 ;; \
    esac; \
    micromamba install -y -n base -c conda-forge -c bioconda "metamdbg=$METAMDBG_VER" \
    && micromamba clean --all --yes

RUN rm -f /opt/conda/lib/python*/site-packages/pip/_vendor/bom.cdx.json \
          /opt/conda/lib/python*/site-packages/pip/_vendor/vendor.txt

ENV PATH="/opt/conda/bin:$PATH"
COPY LICENSE /opt/hvp-lr/LICENSE
# multi-tool image: inherit micromamba base entrypoint
