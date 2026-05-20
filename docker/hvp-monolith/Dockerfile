FROM mambaorg/micromamba:2.4.0-ubuntu24.04

# Each env file is COPY'd immediately before its create step so that editing
# one env file only invalidates that env's layer cache, not the others.

# Main env (base): modern stack on python=3.12 (checkm2, genomad, assemblers,
# binners including concoct, etc.).
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base \
        --channel-priority strict \
        -f /tmp/env.yaml

# Sidecar env (viral): virsorter (>=2) + vibrant. Requires python<=3.10.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-viral.yaml /tmp/env-viral.yaml
RUN micromamba create -y -n viral \
        --channel-priority strict \
        -f /tmp/env-viral.yaml

# Sidecar env (deepvirfinder): python=3.6 + h5py + diamond.
# DeepVirFinder itself is a git repo (jessieren/DeepVirFinder), not a conda
# package; we clone it after the conda solve. theano=1.0.3 + keras=2.2.4 are
# pip-installed because their py3.6 conda builds are not reliable.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-deepvirfinder.yaml /tmp/env-deepvirfinder.yaml
RUN micromamba create -y -n deepvirfinder \
        --channel-priority strict \
        -f /tmp/env-deepvirfinder.yaml && \
    micromamba run -n deepvirfinder pip install --no-cache-dir \
        theano==1.0.3 keras==2.2.4

# Sidecar env (vcontact2): pinned to the working external install referenced
# in the original tpshea2 environment notes.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-vcontact2.yaml /tmp/env-vcontact2.yaml
RUN micromamba create -y -n vcontact2 \
        --channel-priority strict \
        -f /tmp/env-vcontact2.yaml && \
    micromamba clean --all --yes

# DeepVirFinder source. Cloned to /opt/deepvirfinder and exposed via a `dvf.py`
# wrapper in /usr/local/bin (see USER root block below).
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends git ca-certificates && \
    rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 https://github.com/jessieren/DeepVirFinder /opt/deepvirfinder

# Expose every sidecar env's CLIs at the same names as if they were in base
# PATH. Wrapper shims forward args directly to the right env's bin/ so callers
# never need to know which env a tool lives in.
RUN set -e; \
    # viral env wrappers
    for tool in virsorter VIBRANT_run.py VIBRANT_annotation.py VIBRANT_extract_nucleotide.py VIBRANT_extract_protein.py; do \
        printf '#!/bin/bash\nexec /opt/conda/envs/viral/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/$tool"; \
        chmod +x "/usr/local/bin/$tool"; \
    done; \
    # vcontact2 env wrappers
    for tool in vcontact2 vcontact2_gene2genome; do \
        printf '#!/bin/bash\nexec /opt/conda/envs/vcontact2/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/$tool"; \
        chmod +x "/usr/local/bin/$tool"; \
    done; \
    # DeepVirFinder wrapper: activates deepvirfinder env and runs the cloned script.
    printf '#!/bin/bash\nexec micromamba run -n deepvirfinder python /opt/deepvirfinder/dvf.py "$@"\n' \
        > /usr/local/bin/dvf.py; \
    chmod +x /usr/local/bin/dvf.py
# /data is created via WORKDIR below as root by default; make it writable
# for the unprivileged mambauser so tools like run_MaxBin.pl can write
# scratch files in CWD.
RUN mkdir -p /data && chown -R $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER

# Enable conda environment activation for subsequent RUN commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# /usr/local/bin first so the sidecar-env wrappers shadow nothing else.
ENV PATH="/usr/local/bin:/opt/conda/bin:$PATH"

# Verify key tools are available across all envs.
RUN python --version && \
    samtools --version | head -1 && \
    fastqc --version && \
    multiqc --version && \
    concoct --version && \
    virsorter --version && \
    VIBRANT_run.py -h >/dev/null && \
    which vcontact2 && \
    dvf.py -h >/dev/null && \
    test -f /opt/deepvirfinder/dvf.py

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image.
COPY LICENSE /opt/hvp-lr/LICENSE

# Set working directory
WORKDIR /data
