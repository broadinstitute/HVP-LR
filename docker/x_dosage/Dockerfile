# x_dosage — sex-chromosome karyotyping from coverage (T2T-CHM13v2.0).
#
# Bakes the x_dosage scripts (rx_sex.py classifier + calibration/utility tools)
# and the calibrated per-platform configs into a minimal Python image. This image
# is the maintenance home for those scripts (formerly analysis/x_dosage/).
#
# Deps are deliberately minimal: the scripts are stdlib-only (tomllib needs
# Python >=3.11), so there are NO pip/conda installs — just the interpreter,
# the scripts, and the OS utilities the WDL resource-preamble expects.
#
# Multi-arch: linux/amd64 + linux/arm64 (python:3.12-slim is glibc, both arches).
ARG X_DOSAGE_VERSION=0.1.0

FROM --platform=$TARGETPLATFORM python:3.12-slim-bookworm AS runtime
ARG X_DOSAGE_VERSION

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.url="https://github.com/broadinstitute/HVP-LR" \
      org.opencontainers.image.version="${X_DOSAGE_VERSION}" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.description="HVP x_dosage sex-chromosome karyotype caller (T2T-CHM13v2.0)" \
      org.opencontainers.image.vendor="Broad Institute HVP-LR"

# OS utilities: procps for `free -g` in the WDL resource preamble; bash/gawk/
# coreutils/gzip for the classify task's parsing. No compiler or pip needed.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      bash grep gawk coreutils procps gzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Canonical scripts + calibrated configs. This is the maintenance home.
COPY src/ /opt/x_dosage/
RUN chmod +x /opt/x_dosage/*.py
ENV PATH="/opt/x_dosage:${PATH}"

# MIT (repo LICENSE) requires the notice travel with copies. The x_dosage code
# itself is BSD-3-Clause (see image label).
COPY LICENSE /opt/hvp-lr/LICENSE

# No ENTRYPOINT: this image ships several scripts (rx_sex.py, calibrate.py, ...);
# callers invoke the one they need (e.g. `rx_sex.py --config ... in.bed.gz`).
