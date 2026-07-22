#!/usr/bin/env bash
# Pipeline step 12 — sample-structure diagnostic (PCA/UMAP/dendrogram on
# sample×virus) + sample-colored unified viral UMAP.
# Wraps `python -m hvp_viral_viz.cohort_extras`.
#
# Consumes out/<cohort>/_pooled/ produced by step 11. Writes diagnostic
# PNGs + sample_diag_coords.parquet into out/<cohort>/_pooled/plots/.
#
# Usage (from package root):
#   bash scripts/12_cohort_extras.sh --pooled-dir out/cohort_2026-06-22/_pooled all
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.cohort_extras "$@"
