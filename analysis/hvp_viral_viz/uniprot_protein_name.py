"""Fetch UniProt protein_name field for BFVD uniprots.

UniProt's `protein_name` carries the curated recommended/submitted protein
name (e.g. "Capsid protein VP1", "RNA-directed RNA polymerase"). Used to
annotate cluster-marker ORFs with actual protein function, so cluster
labels can be cross-checked against protein content rather than just
host/taxonomy joins.

Outputs:
    refs/cache/uniprot_protein_name.parquet
        uniprot_id : str
        protein_name : str

Cache is reused unless --refresh is given.
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests

from hvp_viral_viz.uniprot_host import (
    UNIPROTKB_RE,
    collect_sample_uniprots,
)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100


def fetch_batch(accessions: list[str], session: requests.Session) -> pd.DataFrame:
    """One batch. Returns uniprot_id, protein_name."""
    query = " OR ".join(f"accession:{a}" for a in accessions)
    params = {
        "format": "tsv",
        "fields": "accession,protein_name",
        "query": query,
        "size": "500",
    }
    r = session.get(UNIPROT_URL, params=params, timeout=60)
    r.raise_for_status()
    lines = r.text.splitlines()
    if len(lines) < 2:
        return pd.DataFrame(columns=["uniprot_id", "protein_name"])
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < 2:
            rows.append((parts[0] if parts else "", ""))
        else:
            rows.append((parts[0], parts[1]))
    return pd.DataFrame(rows, columns=["uniprot_id", "protein_name"])


def fetch_all(accessions: list[str], cache_path: Path | None = None) -> pd.DataFrame:
    session = requests.Session()
    session.headers["User-Agent"] = "hvp-viral-viz/0.1 uniprot-protein-name fetcher"
    results = []
    n = len(accessions)
    t0 = time.time()
    for i in range(0, n, BATCH_SIZE):
        batch = accessions[i : i + BATCH_SIZE]
        try:
            df = fetch_batch(batch, session)
        except requests.HTTPError as e:
            code = e.response.status_code
            if code == 400:
                print(f"[uniprot-pn] batch {i}: 400, bisecting {len(batch)}", file=sys.stderr)
                mid = len(batch) // 2 or 1
                df_a = fetch_batch(batch[:mid], session) if batch[:mid] else pd.DataFrame(columns=["uniprot_id","protein_name"])
                df_b = fetch_batch(batch[mid:], session) if batch[mid:] else pd.DataFrame(columns=["uniprot_id","protein_name"])
                df = pd.concat([df_a, df_b], ignore_index=True)
            else:
                print(f"[uniprot-pn] batch {i}: HTTP {code}, retrying once", file=sys.stderr)
                time.sleep(2)
                df = fetch_batch(batch, session)
        results.append(df)
        if (i // BATCH_SIZE) % 25 == 0:
            elapsed = time.time() - t0
            rate = (i + len(batch)) / max(elapsed, 0.1)
            eta = (n - i - len(batch)) / max(rate, 0.1)
            print(
                f"[uniprot-pn] {i + len(batch):,}/{n:,} ({100 * (i + len(batch)) / n:.1f}%) "
                f"rate={rate:.0f}/s eta={eta:.0f}s",
                file=sys.stderr,
            )
        if (i // BATCH_SIZE) % 50 == 49:
            time.sleep(1)
    out = pd.concat(results, ignore_index=True)
    print(f"[uniprot-pn] done: {len(out):,} rows in {time.time() - t0:.0f}s", file=sys.stderr)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_parquet(cache_path, index=False)
        print(f"[uniprot-pn] cached → {cache_path}", file=sys.stderr)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument("--out-dir", type=Path, default=Path("/workspace/viral-viz/out"))
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--source", choices=["samples", "bfvd"], default="samples")
    args = ap.parse_args()

    cache_path = args.refs_dir / "cache" / "uniprot_protein_name.parquet"

    if args.source == "samples":
        accessions = collect_sample_uniprots(args.out_dir)
    else:
        u = pd.read_parquet(args.refs_dir / "uniprot_taxid.parquet")
        accessions = u["uniprot_id"].dropna().unique().tolist()

    accessions = [a for a in accessions if UNIPROTKB_RE.match(a)]
    print(f"[uniprot-pn] candidate UniProtKB accessions: {len(accessions):,}", file=sys.stderr)

    if cache_path.exists() and not args.refresh:
        cached = pd.read_parquet(cache_path)
        already = set(cached["uniprot_id"])
        accessions = [a for a in accessions if a not in already]
        print(f"[uniprot-pn]   {len(already):,} cached; fetching {len(accessions):,} new", file=sys.stderr)
    else:
        cached = None

    if not accessions:
        print(f"[uniprot-pn] nothing to fetch; cache up-to-date", file=sys.stderr)
        return

    new = fetch_all(accessions, cache_path=None)
    if cached is not None:
        combined = pd.concat([cached, new], ignore_index=True).drop_duplicates(subset=["uniprot_id"], keep="last")
    else:
        combined = new

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(cache_path, index=False)
    n_with_name = combined["protein_name"].astype(bool).sum()
    print(f"[uniprot-pn] cache → {cache_path}  ({len(combined):,} total, {n_with_name:,} with name)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
