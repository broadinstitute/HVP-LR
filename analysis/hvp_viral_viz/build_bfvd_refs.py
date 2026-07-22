"""Collapse raw BFVD TSVs into per-uniprot / per-taxid parquet tables.

Inputs (default /workspace/bfvd/):
    bfvd_metadata.tsv                       — 6 col: uniprot, model_id, avg_plddt, ptm, splitted, version
    bfvd_taxid.tsv                          — 2 col: model_pdb_filename, taxid
    bfvd_taxid_rank_scientificname_lineage.tsv  — 5 col: model_pdb_filename, taxid, rank, scientific_name, lineage

Outputs (refs/):
    uniprot_taxid.parquet      — uniprot_id, taxid
    taxid_lineage.parquet      — taxid, scientific_name, realm, kingdom, phylum, class, order, family, genus, species
    uniprot_qc.parquet         — uniprot_id, model_id, avg_plddt, ptm, splitted, version (best avg_plddt per uniprot when split)
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

LINEAGE_PREFIX = {
    "d_": "realm",
    "k_": "kingdom",
    "p_": "phylum",
    "c_": "class",
    "o_": "order",
    "f_": "family",
    "g_": "genus",
    "s_": "species",
}
RANK_COLS = list(LINEAGE_PREFIX.values())


def parse_lineage(lineage: str) -> dict[str, str | None]:
    """Parse 'd_Viruses;-_Duplodnaviria;k_Heunggongvirae;...' → dict per rank."""
    out: dict[str, str | None] = {r: None for r in RANK_COLS}
    if not isinstance(lineage, str) or not lineage:
        return out
    for tok in lineage.split(";"):
        if len(tok) < 2 or tok[1] != "_":
            continue
        prefix = tok[:2]
        value = tok[2:].strip()
        rank = LINEAGE_PREFIX.get(prefix)
        if rank and value:
            out[rank] = value
    return out


_MODEL_RE = re.compile(r"^(?P<stem>[A-Z0-9]+)(?:_(?P<part>\d+))?")


def parse_model_filename(name: str) -> tuple[str, str]:
    """'A0A2P1GMZ4_2_unrelaxed_rank_001_...pdb' → ('A0A2P1GMZ4_2', 'A0A2P1GMZ4').
    'V9SF81_unrelaxed_...pdb' → ('V9SF81', 'V9SF81').

    Returns (model_id, uniprot_stem).
    """
    stem_unrelaxed = name.split("_unrelaxed", 1)[0]
    m = _MODEL_RE.match(stem_unrelaxed)
    if not m:
        return stem_unrelaxed, stem_unrelaxed
    uniprot = m.group("stem")
    return stem_unrelaxed, uniprot


def build_refs(bfvd_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[refs] reading {bfvd_dir / 'bfvd_taxid_rank_scientificname_lineage.tsv'}", file=sys.stderr)
    lineage_df = pd.read_csv(
        bfvd_dir / "bfvd_taxid_rank_scientificname_lineage.tsv",
        sep="\t",
        header=None,
        names=["model_pdb", "taxid", "rank", "scientific_name", "lineage"],
        dtype={"taxid": "Int64", "rank": "string", "scientific_name": "string", "lineage": "string"},
    )
    print(f"[refs]   {len(lineage_df):,} rows", file=sys.stderr)

    # Per-uniprot taxid table
    print("[refs] parsing model filenames → (model_id, uniprot)", file=sys.stderr)
    parsed = lineage_df["model_pdb"].map(parse_model_filename)
    lineage_df["model_id"] = parsed.map(lambda x: x[0])
    lineage_df["uniprot_id"] = parsed.map(lambda x: x[1])

    uniprot_taxid = (
        lineage_df[["uniprot_id", "taxid"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    n_dup_uni = uniprot_taxid["uniprot_id"].duplicated().sum()
    if n_dup_uni:
        print(f"[refs] WARNING: {n_dup_uni} uniprots have multiple taxids — keeping first", file=sys.stderr)
        uniprot_taxid = uniprot_taxid.drop_duplicates(subset=["uniprot_id"]).reset_index(drop=True)
    uniprot_taxid.to_parquet(out_dir / "uniprot_taxid.parquet", index=False)
    print(f"[refs]   uniprot_taxid.parquet: {len(uniprot_taxid):,} rows", file=sys.stderr)

    # Per-taxid lineage table
    print("[refs] parsing lineage strings → ranks", file=sys.stderr)
    per_taxid = (
        lineage_df[["taxid", "scientific_name", "lineage"]]
        .drop_duplicates(subset=["taxid"])
        .reset_index(drop=True)
    )
    ranks_df = pd.DataFrame.from_records(per_taxid["lineage"].fillna("").map(parse_lineage))
    taxid_lineage = pd.concat(
        [per_taxid[["taxid", "scientific_name"]].reset_index(drop=True), ranks_df],
        axis=1,
    )
    taxid_lineage.to_parquet(out_dir / "taxid_lineage.parquet", index=False)
    coverage = {r: taxid_lineage[r].notna().mean() for r in RANK_COLS}
    print(f"[refs]   taxid_lineage.parquet: {len(taxid_lineage):,} rows", file=sys.stderr)
    print(f"[refs]   rank coverage: " + ", ".join(f"{r}={c:.0%}" for r, c in coverage.items()), file=sys.stderr)

    # Per-uniprot QC table (avg_plddt etc)
    print(f"[refs] reading {bfvd_dir / 'bfvd_metadata.tsv'}", file=sys.stderr)
    meta = pd.read_csv(
        bfvd_dir / "bfvd_metadata.tsv",
        sep="\t",
        header=None,
        names=["uniprot_id", "model_id", "avg_plddt", "ptm", "splitted", "version"],
        dtype={
            "uniprot_id": "string",
            "model_id": "string",
            "avg_plddt": "float64",
            "ptm": "float64",
            "splitted": "Int64",
            "version": "string",
        },
    )
    print(f"[refs]   {len(meta):,} model rows", file=sys.stderr)
    # For split uniprots take best avg_plddt
    uniprot_qc = (
        meta.sort_values("avg_plddt", ascending=False)
        .drop_duplicates(subset=["uniprot_id"], keep="first")
        .reset_index(drop=True)
    )
    uniprot_qc.to_parquet(out_dir / "uniprot_qc.parquet", index=False)
    print(f"[refs]   uniprot_qc.parquet: {len(uniprot_qc):,} unique uniprots", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bfvd-dir", type=Path, default=Path("/workspace/bfvd"))
    ap.add_argument("--out-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    args = ap.parse_args()
    build_refs(args.bfvd_dir, args.out_dir)


if __name__ == "__main__":
    main()
