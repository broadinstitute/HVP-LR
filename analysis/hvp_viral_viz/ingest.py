"""Ingest foldseek m8 hit tables into AnnData on disk.

Single-sample mode:  ingest one m8, produce per-virus hit summary + AnnData
                     where obs is one row (the sample) and var is one row per
                     virus (taxid by default).

Cohort mode (future): ingest many m8 files, stack into a samples × viruses
                      count matrix.  Stub for now — same builder functions.

Schema (m8 — foldseek default 7-col output, no header):
    0 query    str
    1 target   str   (BFVD uniprot or split-stem like A0A2P1GMZ4_2)
    2 evalue   float
    3 bits     int
    4 fident   float
    5 alnlen   int
    6 mismatch int

Outputs (one per sample):
    out/<sample>/anndata.h5ad
    out/<sample>/hits_filtered.parquet
    out/<sample>/summary.json
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Discovery / high-confidence thresholds (see THRESHOLD.md)
DISCOVERY = {"evalue_max": 1e-5, "bits_min": 50, "alnlen_min": 50}
HIGHCONF = {"evalue_max": 1e-10, "bits_min": 300, "alnlen_min": 80}

M8_COLS = ["query", "target", "evalue", "bits", "fident", "alnlen", "mismatch"]
M8_DTYPES = {
    "query": "string",
    "target": "string",
    "evalue": "float64",
    "bits": "Int64",
    "fident": "float64",
    "alnlen": "Int64",
    "mismatch": "Int64",
}


def _strip_split_suffix(target: str) -> str:
    """'A0A2P1GMZ4_2' (BFVD split-domain) → 'A0A2P1GMZ4'.

    BFVD split targets carry _N suffix where N is the split index. The
    uniprot accession is everything before the first _N. Pure UniProt
    accessions do not contain underscores.
    """
    if not isinstance(target, str):
        return target
    if "_" not in target:
        return target
    head, tail = target.rsplit("_", 1)
    if tail.isdigit():
        return head
    return target


def passes(df: pd.DataFrame, thr: dict) -> pd.Series:
    return (
        (df["evalue"] <= thr["evalue_max"])
        & (df["bits"] >= thr["bits_min"])
        & (df["alnlen"] >= thr["alnlen_min"])
    )


def parse_orf_source(query: str) -> str:
    """Tag ORF source from query name prefix.

    The HVP pipeline concatenates 3 ORF callers' output, each prefixed
    in the query name: 'vs2|...', 'assembly|...', 'rescued|...'.
    """
    if not isinstance(query, str):
        return "unknown"
    head = query.split("|", 1)[0]
    return head if head in ("vs2", "assembly", "rescued", "genomad") else "unknown"


def load_refs(refs_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load BFVD/host reference tables built by build_bfvd_refs + build_host_table."""
    uniprot_taxid = pd.read_parquet(refs_dir / "uniprot_taxid.parquet")
    taxid_lineage = pd.read_parquet(refs_dir / "taxid_lineage.parquet")
    uniprot_qc = pd.read_parquet(refs_dir / "uniprot_qc.parquet")
    taxid_host = pd.read_parquet(refs_dir / "taxid_host.parquet")
    return uniprot_taxid, taxid_lineage, uniprot_qc, taxid_host


def annotate_hits(
    m8: pd.DataFrame,
    uniprot_taxid: pd.DataFrame,
    taxid_lineage: pd.DataFrame,
    uniprot_qc: pd.DataFrame,
    taxid_host: pd.DataFrame,
) -> pd.DataFrame:
    """Join hits to BFVD metadata, lineage, host."""
    df = m8.copy()
    df["uniprot_id"] = df["target"].map(_strip_split_suffix).astype("string")
    df["orf_source"] = df["query"].map(parse_orf_source).astype("string")

    df = df.merge(uniprot_taxid, on="uniprot_id", how="left")
    df = df.merge(taxid_lineage, on="taxid", how="left")
    df = df.merge(
        uniprot_qc[["uniprot_id", "avg_plddt", "ptm", "splitted", "version"]],
        on="uniprot_id",
        how="left",
    )
    df = df.merge(taxid_host, on="taxid", how="left")
    return df


def build_count_matrix(
    hits: pd.DataFrame,
    sample: str,
    grouping: str = "taxid",
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Build samples × viruses count matrix.

    Each passing (query, target) row contributes one count to its virus.
    `grouping` ∈ {'taxid', 'uniprot', 'species', 'genus', 'family'}.

    Returns (X, obs, var) where:
        X.shape = (n_samples=1, n_viruses)
        obs.index = [sample]
        var.index = grouping values (taxids as strings)
    """
    if grouping == "uniprot":
        grp_col = "uniprot_id"
    else:
        grp_col = grouping  # 'taxid' / 'species' / 'genus' / 'family'

    counts = (
        hits.groupby(grp_col, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    # Build var with metadata; one row per virus group
    if grouping == "taxid":
        meta_cols = ["taxid", "scientific_name", "realm", "kingdom", "phylum", "class",
                     "order", "family", "genus", "species",
                     "host_group", "host_source", "host_confidence"]
        meta = (
            hits[meta_cols]
            .drop_duplicates(subset=["taxid"])
            .reset_index(drop=True)
        )
        var = counts.merge(meta, on="taxid", how="left")
        var["taxid_str"] = var["taxid"].astype(str)
        var = var.set_index("taxid_str")
    elif grouping == "uniprot":
        meta_cols = ["uniprot_id", "taxid", "scientific_name", "family", "genus",
                     "species", "host_group", "host_confidence",
                     "avg_plddt", "ptm", "splitted", "version"]
        meta = hits[meta_cols].drop_duplicates(subset=["uniprot_id"]).reset_index(drop=True)
        var = counts.merge(meta, on="uniprot_id", how="left").set_index("uniprot_id")
    else:
        var = counts.set_index(grp_col)

    # Build obs (single row for single-sample mode)
    obs = pd.DataFrame({"sample": [sample]}).set_index("sample")

    X = var["count"].astype(np.int32).values.reshape(1, -1)
    return X, obs, var


def ingest_single(
    m8_path: Path,
    sample: str,
    refs_dir: Path,
    out_dir: Path,
    grouping: str = "taxid",
) -> None:
    """Ingest one sample's m8 into AnnData on disk."""
    import anndata as ad

    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"[ingest] {sample}: reading {m8_path}", file=sys.stderr)
    m8 = pd.read_csv(m8_path, sep="\t", header=None, names=M8_COLS, dtype=M8_DTYPES)
    print(f"[ingest]   {len(m8):,} raw hits, {time.time() - t0:.1f}s", file=sys.stderr)

    # Apply discovery threshold
    mask_d = passes(m8, DISCOVERY)
    mask_hc = passes(m8, HIGHCONF)
    n_total = len(m8)
    n_disc = int(mask_d.sum())
    n_hc = int(mask_hc.sum())
    print(
        f"[ingest]   discovery passing: {n_disc:,} ({100 * n_disc / n_total:.1f}%); "
        f"high-conf: {n_hc:,} ({100 * n_hc / n_total:.1f}%)",
        file=sys.stderr,
    )

    hits = m8.loc[mask_d].copy()
    hits["is_high_conf"] = mask_hc.loc[mask_d].values

    # Load refs + annotate
    uniprot_taxid, taxid_lineage, uniprot_qc, taxid_host = load_refs(refs_dir)
    hits = annotate_hits(hits, uniprot_taxid, taxid_lineage, uniprot_qc, taxid_host)

    # Join coverage diagnostics
    n_unmapped_uniprot = hits["taxid"].isna().sum()
    print(
        f"[ingest]   unmapped uniprot→taxid: {n_unmapped_uniprot:,} "
        f"({100 * n_unmapped_uniprot / len(hits):.2f}%)",
        file=sys.stderr,
    )

    # Persist annotated hits
    hits.to_parquet(out_dir / "hits_filtered.parquet", index=False)
    print(f"[ingest]   wrote hits_filtered.parquet: {len(hits):,} rows", file=sys.stderr)

    # Build count matrix and AnnData
    X, obs, var = build_count_matrix(hits, sample=sample, grouping=grouping)
    print(f"[ingest]   count matrix: {X.shape[0]}×{X.shape[1]} ({grouping})", file=sys.stderr)

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Stash threshold + summary as uns
    adata.uns["thresholds"] = {"discovery": DISCOVERY, "high_conf": HIGHCONF}
    adata.uns["counts"] = {
        "n_raw": int(n_total),
        "n_discovery": int(n_disc),
        "n_high_conf": int(n_hc),
        "n_unmapped_uniprot": int(n_unmapped_uniprot),
    }
    adata.uns["grouping"] = grouping

    out_h5ad = out_dir / "anndata.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"[ingest]   wrote {out_h5ad} ({out_h5ad.stat().st_size / 1e6:.1f} MB)", file=sys.stderr)

    # Summary JSON
    by_host = (
        var.assign(c=var["count"])
        .groupby("host_group", dropna=False)["c"]
        .sum()
        .to_dict()
    )
    by_realm = (
        var.assign(c=var["count"])
        .groupby("realm", dropna=False)["c"]
        .sum()
        .to_dict()
    )
    by_orf_source = hits["orf_source"].value_counts().to_dict()

    summary = {
        "sample": sample,
        "input_m8": str(m8_path),
        "thresholds": {"discovery": DISCOVERY, "high_conf": HIGHCONF},
        "counts": {
            "n_raw_hits": int(n_total),
            "n_discovery_hits": int(n_disc),
            "n_high_conf_hits": int(n_hc),
            "n_unique_queries_with_hit": int(hits["query"].nunique()),
            "n_unique_targets": int(hits["target"].nunique()),
            "n_unique_uniprots": int(hits["uniprot_id"].nunique()),
            "n_unique_taxids": int(hits["taxid"].nunique()),
            "n_unique_families": int(hits["family"].nunique()),
        },
        "hits_by_host_group": {str(k): int(v) for k, v in by_host.items()},
        "hits_by_realm": {str(k): int(v) for k, v in by_realm.items()},
        "hits_by_orf_source": {str(k): int(v) for k, v in by_orf_source.items()},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"[ingest]   wrote summary.json", file=sys.stderr)
    print(f"[ingest] done in {time.time() - t0:.1f}s", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--m8", type=Path, required=True, help="foldseek m8 hits TSV")
    ap.add_argument("--sample", required=True, help="sample id (used as obs key + out subdir)")
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument("--out-dir", type=Path, default=Path("/workspace/viral-viz/out"))
    ap.add_argument("--grouping", default="taxid",
                    choices=["taxid", "uniprot", "species", "genus", "family"])
    args = ap.parse_args()

    out_dir = args.out_dir / args.sample
    ingest_single(args.m8, args.sample, args.refs_dir, out_dir, args.grouping)


if __name__ == "__main__":
    main()
