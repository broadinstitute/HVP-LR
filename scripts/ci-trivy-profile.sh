#!/usr/bin/env bash
# Pick a Trivy scan profile based on the local image's uncompressed size.
#
# Big scientific images (multi-env conda + model weights + reference DBs) trip
# the secret scanner's per-file size guard and blow past Trivy's default 5-min
# scan deadline (observed: "context deadline exceeded" on a ~5 GB image). The
# secret scanner is designed for source repositories with API keys, not for
# scientific containers full of binary model weights — dropping it on large
# images is safe.
#
# Smaller images get the full default profile (vuln + secret) at the default
# timeout, so generic helper images aren't downgraded just because they share
# this workflow.
#
# Driven by:
#   IMAGE_REF                  : image ref already loaded locally (required)
#   TRIVY_LARGE_IMAGE_BYTES    : threshold in bytes (default 2 GiB)
#
# Emits to $GITHUB_OUTPUT:
#   scanners=vuln              | vuln,secret
#   timeout=20m                | 10m
#   size_bytes=<int>           (for logs)
set -euo pipefail

: "${IMAGE_REF:?IMAGE_REF must be set}"
threshold="${TRIVY_LARGE_IMAGE_BYTES:-2147483648}"

size="$(docker image inspect -f '{{.Size}}' "$IMAGE_REF")"

if [ "$size" -gt "$threshold" ]; then
  scanners="vuln"
  timeout="20m"
  echo "::notice title=Trivy profile::Image $IMAGE_REF is $size bytes (> $threshold); vuln-only scan, 20m timeout."
else
  scanners="vuln,secret"
  timeout="10m"
  echo "::notice title=Trivy profile::Image $IMAGE_REF is $size bytes (<= $threshold); full scan (vuln+secret), 10m timeout."
fi

{
  echo "scanners=$scanners"
  echo "timeout=$timeout"
  echo "size_bytes=$size"
} >> "$GITHUB_OUTPUT"
