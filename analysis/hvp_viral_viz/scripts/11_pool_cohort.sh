#!/usr/bin/env bash
# Pipeline step 11 — pool all per-sample hits into a unified virus×pooled-ORF
# landscape UMAP + leiden clusters.
# Wraps `python -m hvp_viral_viz.pool_cohort`.
#
# Consumes out/<cohort>/<sample>/ produced by step 10. Writes
# out/<cohort>/_pooled/{pooled_hits.parquet, pooled_anndata.h5ad,
# virus_samples.parquet, summary.json, plots/}.
#
# Usage (from package root):
#   bash scripts/11_pool_cohort.sh --cohort-dir out/cohort_2026-06-22 [--leiden-resolution N]
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"
exec python -m hvp_viral_viz.pool_cohort "$@"
