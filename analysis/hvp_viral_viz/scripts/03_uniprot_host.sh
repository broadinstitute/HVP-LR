#!/usr/bin/env bash
# Pipeline step 03 — fetch UniProt virus_hosts cross-reference for BFVD uniprots.
# Wraps `python -m hvp_viral_viz.uniprot_host`.
#
# Produces refs/cache/uniprot_host.parquet; cache is reused unless --refresh.
#
# Usage (from package root):
#   bash scripts/03_uniprot_host.sh --source bfvd      # all BFVD uniprots
#   bash scripts/03_uniprot_host.sh --source samples   # only uniprots seen in out/<sample>/
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.uniprot_host "$@"
