#!/usr/bin/env bash
# Pipeline step 05 — build per-taxid host annotation table (4-layer fill).
# Wraps `python -m hvp_viral_viz.build_host_table`.
#
# Consumes refs/ (built by 02), refs/cache/uniprot_host.parquet (03),
# refs/cache/host_taxid_to_group.parquet (04), and the ICTV VMR xlsx
# downloaded by 01. Produces refs/taxid_host.parquet + audit TSV.
#
# Usage (from package root):
#   bash scripts/05_build_host_table.sh [--refs-dir DIR] [--ictv-vmr PATH]
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.build_host_table "$@"
