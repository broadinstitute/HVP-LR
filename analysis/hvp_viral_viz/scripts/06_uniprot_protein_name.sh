#!/usr/bin/env bash
# Pipeline step 06 — fetch UniProt protein_name field for BFVD uniprots.
# Wraps `python -m hvp_viral_viz.uniprot_protein_name`.
#
# Optional but required for label_clusters --with-protein-markers (invoked
# by 10_run_cohort.sh by default). Produces refs/cache/uniprot_protein_name.parquet.
#
# Usage (from package root):
#   bash scripts/06_uniprot_protein_name.sh --source bfvd
#   bash scripts/06_uniprot_protein_name.sh --source samples
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.uniprot_protein_name "$@"
