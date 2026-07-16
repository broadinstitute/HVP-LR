#!/usr/bin/env bash
# Per-sample ingest + plot pipeline for the 45-sample cohort.
#
# For each m8 in data/cohort_2026-06-22/raw/ run:
#   1. ingest        — filter + annotate hits → anndata.h5ad + hits_filtered.parquet
#   2. plots         — tier-1 plots (abundance bar/family, heatmap, UMAPs)
#   3. plots_tier2   — tier-2 plots (sankey/sunburst/bits-fident/plDDT/coverage)
#   4. label_clusters — biological labels per leiden cluster + labeled UMAP
#
# Outputs land in out/cohort_2026-06-22/<sample>/.
# Re-runs are idempotent: per-stage marker file gates re-execution. Delete the
# sample's out subdir to force rebuild.
#
# Usage (from package root):
#   bash scripts/10_run_cohort.sh               # sequential
#   PARALLEL=4 bash scripts/10_run_cohort.sh    # 4 samples in parallel

set -euo pipefail

# Land in analysis/hvp_viral_viz/ (script's parent) so paths below are relative.
# Set PYTHONPATH to its parent so `python -m hvp_viral_viz.X` resolves the package.
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)/..${PYTHONPATH:+:$PYTHONPATH}"

RAW=data/cohort_2026-06-22/raw
OUT=out/cohort_2026-06-22
REFS=refs
LOGS=$OUT/_logs
PARALLEL=${PARALLEL:-1}

mkdir -p "$OUT" "$LOGS"

run_one() {
    local m8="$1"
    local sample
    sample=$(basename "$m8" .vs_bfvd.m8)
    local sample_dir="$OUT/$sample"
    local log="$LOGS/$sample.log"

    mkdir -p "$sample_dir"
    {
        echo "=== $(date -Is) $sample begin ==="

        if [[ ! -f "$sample_dir/anndata.h5ad" ]]; then
            echo "--- ingest ---"
            python -m hvp_viral_viz.ingest \
                --m8 "$m8" --sample "$sample" \
                --refs-dir "$REFS" --out-dir "$OUT"
        else
            echo "--- ingest [skip: anndata.h5ad exists] ---"
        fi

        if [[ ! -f "$sample_dir/plots/umap_by_family.png" ]]; then
            echo "--- plots (tier 1) ---"
            python -m hvp_viral_viz.plots --sample-dir "$sample_dir"
        else
            echo "--- plots [skip: umap_by_family.png exists] ---"
        fi

        if [[ ! -f "$sample_dir/plots/sankey_orf_to_host.html" ]]; then
            echo "--- plots_tier2 ---"
            python -m hvp_viral_viz.plots_tier2 --sample-dir "$sample_dir"
        else
            echo "--- plots_tier2 [skip: sankey_orf_to_host.html exists] ---"
        fi

        if [[ ! -f "$sample_dir/plots/umap_by_cluster_labeled_order_with_proteins.png" ]]; then
            echo "--- label_clusters ---"
            python -m hvp_viral_viz.label_clusters \
                --sample-dir "$sample_dir" --rank order \
                --refs-dir "$REFS" --with-protein-markers
        else
            echo "--- label_clusters [skip: labeled UMAP exists] ---"
        fi

        echo "=== $(date -Is) $sample done ==="
    } &> "$log"
    echo "[run_cohort] $sample → $log"
}

export -f run_one
export OUT REFS LOGS

if [[ "$PARALLEL" -gt 1 ]]; then
    echo "[run_cohort] parallel=$PARALLEL"
    find "$RAW" -name '*.vs_bfvd.m8' | sort | xargs -n1 -P"$PARALLEL" -I{} bash -c 'run_one "$@"' _ {}
else
    echo "[run_cohort] sequential"
    for m8 in $(find "$RAW" -name '*.vs_bfvd.m8' | sort); do
        run_one "$m8"
    done
fi

echo "[run_cohort] all done"
