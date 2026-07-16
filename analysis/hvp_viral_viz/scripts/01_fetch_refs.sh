#!/usr/bin/env bash
# Fetch the backing data sources required by build_bfvd_refs.py + build_host_table.py.
#
# Downloads (idempotent — skips files already present and non-empty):
#   1. BFVD metadata + taxid + lineage TSVs from steineggerlab.
#   2. ICTV VMR MSL41 xlsx from ictv.global.
#
# Pipeline step 01. Run from the package root (analysis/hvp_viral_viz/).
#
# Usage:
#   bash scripts/01_fetch_refs.sh [--out-dir DIR]
#
# Next steps run by sibling 02_/05_ wrappers — see README.md "Pipeline".

set -euo pipefail

OUT_DIR="./refs/downloads"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$OUT_DIR/bfvd"

# ---- BFVD ----------------------------------------------------------------
BFVD_BASE="https://bfvd.steineggerlab.workers.dev/latest"
BFVD_FILES=(
    "bfvd_metadata.tsv"
    "bfvd_taxid.tsv"
    "bfvd_taxid_rank_scientificname_lineage.tsv"
)

for f in "${BFVD_FILES[@]}"; do
    dest="$OUT_DIR/bfvd/$f"
    if [[ -s "$dest" ]]; then
        echo "[fetch_refs] skip BFVD $f (already present)"
        continue
    fi
    echo "[fetch_refs] downloading BFVD $f"
    wget -c -O "$dest.partial" "$BFVD_BASE/$f"
    mv "$dest.partial" "$dest"
done

# ---- ICTV VMR ------------------------------------------------------------
# Canonical filename build_host_table.py looks for: ictv_vmr_msl41.xlsx.
# Upstream filename embeds release date and may change between MSLs; we
# rename to the stable local name.
ICTV_URL="https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260320.xlsx"
ICTV_DEST="$OUT_DIR/ictv_vmr_msl41.xlsx"

if [[ -s "$ICTV_DEST" ]]; then
    echo "[fetch_refs] skip ICTV VMR (already present)"
else
    echo "[fetch_refs] downloading ICTV VMR MSL41"
    wget -c -O "$ICTV_DEST.partial" "$ICTV_URL"
    mv "$ICTV_DEST.partial" "$ICTV_DEST"
fi

# ---- Summary -------------------------------------------------------------
echo
echo "[fetch_refs] done. files in $OUT_DIR/:"
ls -lh "$OUT_DIR/bfvd/" "$OUT_DIR/"*.xlsx 2>/dev/null || true
echo
echo "Next steps:"
echo "  python -m hvp_viral_viz.build_bfvd_refs  --bfvd-dir $OUT_DIR/bfvd  --out-dir refs"
echo "  python -m hvp_viral_viz.build_host_table --refs-dir refs --ictv-vmr $ICTV_DEST"
