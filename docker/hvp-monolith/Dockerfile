# ----------------------------------------------------------------------------
# Stage 1 (builder): assemble every conda env, clone DeepVirFinder, write the
# /usr/local/bin wrapper shims, then aggressively strip caches and docs from
# /opt/conda before the final stage copies it.
# ----------------------------------------------------------------------------
FROM mambaorg/micromamba:2.4.0-ubuntu24.04 AS builder

# Main env (base): modern stack on python=3.12 (checkm2, genomad, assemblers,
# binners including concoct, etc.).
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base \
        --channel-priority strict \
        -f /tmp/env.yaml && \
    micromamba clean --all --yes

# metamdbg (per-arch). bioconda ships metamdbg 1.4 on linux-64 only; the
# latest linux-aarch64 build is 1.2. Pick the newest available per arch
# rather than downgrading amd64 just to keep one yaml line.
ARG TARGETARCH
RUN set -e; \
    case "$TARGETARCH" in \
        amd64) METAMDBG_VER=1.4 ;; \
        arm64) METAMDBG_VER=1.2 ;; \
        *) echo "unsupported TARGETARCH=$TARGETARCH" >&2; exit 1 ;; \
    esac; \
    micromamba install -y -n base \
        --channel-priority strict \
        -c conda-forge -c bioconda \
        "metamdbg=$METAMDBG_VER" && \
    micromamba clean --all --yes

# Sidecar env (viral): virsorter (>=2) + vibrant. Requires python<=3.10.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-viral.yaml /tmp/env-viral.yaml
RUN micromamba create -y -n viral \
        --channel-priority strict \
        -f /tmp/env-viral.yaml && \
    micromamba clean --all --yes

# Sidecar env (deepvirfinder): python=3.6 + h5py + biopython + diamond.
# DeepVirFinder itself is a git repo (jessieren/DeepVirFinder), not a conda
# package; we clone it after the conda solve. theano=1.0.3 + keras=2.2.4 are
# pip-installed because their py3.6 conda builds are not reliable.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-deepvirfinder.yaml /tmp/env-deepvirfinder.yaml
RUN micromamba create -y -n deepvirfinder \
        --channel-priority strict \
        -f /tmp/env-deepvirfinder.yaml && \
    micromamba run -n deepvirfinder pip install --no-cache-dir \
        theano==1.0.3 keras==2.2.4 && \
    micromamba clean --all --yes

# Sidecar env (vcontact2): pinned to the working external install referenced
# in the original tpshea2 environment notes.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-vcontact2.yaml /tmp/env-vcontact2.yaml
RUN micromamba create -y -n vcontact2 \
        --channel-priority strict \
        -f /tmp/env-vcontact2.yaml && \
    micromamba clean --all --yes

# Sidecar env (checkm2): tensorflow=2.17 (pip wheel — conda-forge ships no
# linux-aarch64 build at that version) + checkm2 1.1.0 from source.
COPY --chown=$MAMBA_USER:$MAMBA_USER env-checkm2.yaml /tmp/env-checkm2.yaml
RUN micromamba create -y -n checkm2 \
        --channel-priority strict \
        -f /tmp/env-checkm2.yaml && \
    micromamba clean --all --yes

# DeepVirFinder source. Cloned to /opt/deepvirfinder; .git/ stripped because
# the runtime doesn't need git history.
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends git ca-certificates && \
    rm -rf /var/lib/apt/lists/* && \
    git clone --depth 1 https://github.com/jessieren/DeepVirFinder /opt/deepvirfinder && \
    rm -rf /opt/deepvirfinder/.git

# Wrapper shims so callers never have to know which env owns a tool.
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
    # checkm2 env wrapper
    printf '#!/bin/bash\nexec /opt/conda/envs/checkm2/bin/checkm2 "$@"\n' \
        > /usr/local/bin/checkm2; \
    chmod +x /usr/local/bin/checkm2; \
    # DeepVirFinder wrapper
    printf '#!/bin/bash\nexec micromamba run -n deepvirfinder python /opt/deepvirfinder/dvf.py "$@"\n' \
        > /usr/local/bin/dvf.py; \
    chmod +x /usr/local/bin/dvf.py

# Strip caches, byte-compiled modules, and docs from /opt/conda. Python will
# regenerate __pycache__ at runtime for modules it imports; nothing here is
# required at runtime. The /opt/conda/pkgs/* extracted dirs hold a second copy
# of every installed package (~9 GB) plus stale versions left behind by env
# solves — those stale copies trip Trivy with CVEs against versions that are
# not actually on $PATH in any env.
RUN set -e; \
    find /opt/conda -type d -name __pycache__ -prune -exec rm -rf {} +; \
    find /opt/conda -type f -name '*.pyc' -delete; \
    find /opt/conda -type f -name '*.pyo' -delete; \
    find /opt/conda -type d -path '*/share/man*' -prune -exec rm -rf {} +; \
    find /opt/conda -type d -path '*/share/doc/*' -prune -exec rm -rf {} +; \
    find /opt/conda -type d -path '*/share/info/*' -prune -exec rm -rf {} +; \
    rm -rf /opt/conda/pkgs 2>/dev/null || true; \
    rm -rf /home/$MAMBA_USER/.mamba/pkgs 2>/dev/null || true; \
    true

# ----------------------------------------------------------------------------
# Stage 2 (final): minimal mambaorg base + copied conda tree + DVF + wrappers.
# Drops the apt git/ca-certificates layer and the builder's intermediate
# pkg caches that would otherwise live in earlier image layers.
# ----------------------------------------------------------------------------
FROM mambaorg/micromamba:2.4.0-ubuntu24.04

USER root
# COPY --from preserves ownership from the builder stage (already MAMBA_USER
# because every file was created in the builder as that user). Re-issuing
# --chown here is correct but multiplies the COPY time by the file count
# (millions for a conda tree).
COPY --from=builder /opt/conda         /opt/conda
COPY --from=builder /opt/deepvirfinder /opt/deepvirfinder
COPY --from=builder /usr/local/bin     /usr/local/bin
RUN mkdir -p /data && chown -R $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER

# Enable conda environment activation for subsequent RUN commands.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# /usr/local/bin first so the sidecar-env wrappers shadow nothing else.
ENV PATH="/usr/local/bin:/opt/conda/bin:$PATH"

# Smoke test: every env answers at least to --help / --version.
RUN python --version && \
    samtools --version | head -1 && \
    fastqc --version && \
    multiqc --version && \
    concoct --version && \
    virsorter --version && \
    VIBRANT_run.py -h >/dev/null && \
    which vcontact2 && \
    checkm2 --version && \
    dvf.py -h >/dev/null && \
    test -f /opt/deepvirfinder/dvf.py

# MIT requires the copyright notice + license text travel with copies of
# the software. Bundle the repo-root LICENSE into the image.
COPY LICENSE /opt/hvp-lr/LICENSE

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.description="HVP-LR monolith: VirSorter2, geNomad, CheckM2, vConTACT2, DeepVirFinder, VIBRANT, MultiQC, FastQC, samtools, CONCOCT, MetaMDBG, and supporting bioinformatics utilities for viral metagenomics" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

WORKDIR /data
