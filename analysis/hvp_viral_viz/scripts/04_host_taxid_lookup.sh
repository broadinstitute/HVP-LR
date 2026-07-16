#!/usr/bin/env bash
# Pipeline step 04 — map host_taxid → host_group via NCBI eutils.
# Wraps `python -m hvp_viral_viz.host_taxid_lookup`.
#
# Consumes uniprot_host.parquet (produced by step 03), produces
# refs/cache/host_taxid_to_group.parquet (consumed by step 05).
#
# Usage (from package root):
#   bash scripts/04_host_taxid_lookup.sh [--refs-dir DIR]
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.host_taxid_lookup "$@"
