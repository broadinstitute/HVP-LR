# vConTACT2 — viral genome clustering (py3.8 sidecar). Pins verbatim from the
# monolith. py3.8 cannot take setuptools>=78, so its old-setuptools/py CVEs are
# handled via .trivyignore (VEX: no WDL task invokes this tool).
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

USER root
RUN apt-get update && apt-get upgrade -y && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://bitbucket.org/MAVERICLab/vcontact2" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.description="vConTACT2 viral genome clustering (py3.8)" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete

RUN rm -f /opt/conda/lib/python*/site-packages/pip/_vendor/bom.cdx.json \
          /opt/conda/lib/python*/site-packages/pip/_vendor/vendor.txt 2>/dev/null || true

ENV PATH="/opt/conda/bin:$PATH"
COPY LICENSE /opt/hvp-lr/LICENSE
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "vcontact2"]
