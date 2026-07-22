#!/usr/bin/env bash
# Pipeline step 13 — recolor pooled-cohort plots by specimen type (stool vs saliva).
# Wraps `python -m hvp_viral_viz.cohort_specimen`.
#
# Consumes out/<cohort>/_pooled/ (step 11) + sample_diag_coords.parquet
# (step 12) + a Terra data table TSV produced by step 07. Writes
# *_by_specimen.png + umap_per_specimen.png + umap_by_specimen_fraction.png.
#
# Usage (from package root):
#   bash scripts/13_cohort_specimen.sh \
#     --pooled-dir out/cohort_2026-06-22/_pooled \
#     --table-tsv data/terra_tables/2026-06-24/HVP-0006_1.tsv
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.cohort_specimen "$@"
