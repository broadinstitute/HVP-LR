FROM mambaorg/micromamba:2.4.0-ubuntu24.04

# Copy environment file and install dependencies
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Enable conda environment activation for subsequent RUN commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Add manually to PATH for all the image-users who override our entrypoint
ENV PATH="/opt/conda/bin:$PATH"

# Verify key tools are available
RUN python --version && \
    samtools --version | head -1 && \
    fastqc --version && \
    multiqc --version

# Set working directory
WORKDIR /data
