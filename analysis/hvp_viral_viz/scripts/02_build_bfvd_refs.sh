#!/usr/bin/env bash
# Pipeline step 02 — collapse raw BFVD TSVs into per-uniprot / per-taxid parquet tables.
# Wraps `python -m hvp_viral_viz.build_bfvd_refs`.
#
# Run after 01_fetch_refs.sh has populated refs/downloads/bfvd/.
#
# Usage (from package root):
#   bash scripts/02_build_bfvd_refs.sh [--bfvd-dir DIR] [--out-dir DIR]
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.build_bfvd_refs "$@"
