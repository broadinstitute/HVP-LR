# DeepVirFinder — git-cloned tool (jessieren/DeepVirFinder) on a py3.6 sidecar.
# theano/keras are pip-installed (py3.6 conda builds unreliable). Recipe mirrors
# the retired monolith. py3.6 is EOL; its fixed-CVE debt is handled via
# .trivyignore (VEX: no WDL task invokes this tool).
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

USER root
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/jessieren/DeepVirFinder" \
      org.opencontainers.image.licenses="USC-RL v1.0" \
      org.opencontainers.image.description="DeepVirFinder viral sequence identification (py3.6 + theano/keras)" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba run -n base pip install --no-cache-dir theano==1.0.3 keras==2.2.4 \
    && micromamba clean --all --yes \
    && find /opt/conda -follow -type f -name '*.pyc' -delete

# DeepVirFinder source (strip .git — runtime doesn't need history).
USER root
RUN git clone --depth 1 https://github.com/jessieren/DeepVirFinder /opt/deepvirfinder \
    && rm -rf /opt/deepvirfinder/.git \
    && chown -R $MAMBA_USER:$MAMBA_USER /opt/deepvirfinder
# Wrapper so callers just run `dvf.py`.
RUN printf '#!/bin/bash\nexec /opt/conda/bin/python /opt/deepvirfinder/dvf.py "$@"\n' > /usr/local/bin/dvf.py \
    && chmod +x /usr/local/bin/dvf.py
USER $MAMBA_USER

RUN rm -f /opt/conda/lib/python*/site-packages/pip/_vendor/bom.cdx.json \
          /opt/conda/lib/python*/site-packages/pip/_vendor/vendor.txt 2>/dev/null || true

ENV KERAS_BACKEND=theano
ENV PATH="/usr/local/bin:/opt/conda/bin:$PATH"
COPY LICENSE /opt/hvp-lr/LICENSE
ENTRYPOINT ["/usr/local/bin/dvf.py"]
