"""Fetch UniProt virus_hosts field for BFVD uniprots.

UniProt's `virus_hosts` cross-reference carries per-entry host species
(curated for characterized viruses). Format returned by REST as TSV:
  "Homo sapiens (Human) [TaxID: 9606]; Mus musculus (Mouse) [TaxID: 10090]"

Outputs:
    refs/cache/uniprot_host.parquet
        uniprot_id : str
        host_organism : str   (raw string; ';'-joined if multiple)
        host_taxid : str      (first parsed taxid; ';'-joined if multiple)

Cache is reused unless --refresh is given.
"""
from __future__ import annotations
import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100
TAXID_RE = re.compile(r"\[TaxID:\s*(\d+)\]")
# Standard UniProtKB accession pattern (excludes UPI* UniParc IDs)
UNIPROTKB_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def fetch_batch(accessions: list[str], session: requests.Session) -> pd.DataFrame:
    """One batch of up to BATCH_SIZE accessions. Returns Entry, Virus hosts tsv."""
    query = " OR ".join(f"accession:{a}" for a in accessions)
    params = {
        "format": "tsv",
        "fields": "accession,virus_hosts",
        "query": query,
        "size": "500",
    }
    r = session.get(UNIPROT_URL, params=params, timeout=60)
    r.raise_for_status()
    lines = r.text.splitlines()
    if len(lines) < 2:
        return pd.DataFrame(columns=["uniprot_id", "host_organism"])
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < 2:
            rows.append((parts[0] if parts else "", ""))
        else:
            rows.append((parts[0], parts[1]))
    return pd.DataFrame(rows, columns=["uniprot_id", "host_organism"])


def parse_taxids(host_str: str) -> str:
    if not isinstance(host_str, str) or not host_str:
        return ""
    matches = TAXID_RE.findall(host_str)
    return ";".join(matches)


def fetch_all(accessions: list[str], cache_path: Path | None = None) -> pd.DataFrame:
    session = requests.Session()
    session.headers["User-Agent"] = "hvp-viral-viz/0.1 uniprot-host fetcher"
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
                # Bisect — a poisoned accession in this batch. Halve and retry.
                print(f"[uniprot-host] batch {i}: 400, bisecting {len(batch)}", file=sys.stderr)
                mid = len(batch) // 2 or 1
                df_a = fetch_batch(batch[:mid], session) if batch[:mid] else pd.DataFrame(columns=["uniprot_id","host_organism"])
                df_b = fetch_batch(batch[mid:], session) if batch[mid:] else pd.DataFrame(columns=["uniprot_id","host_organism"])
                df = pd.concat([df_a, df_b], ignore_index=True)
            else:
                print(f"[uniprot-host] batch {i}: HTTP {code}, retrying once", file=sys.stderr)
                time.sleep(2)
                df = fetch_batch(batch, session)
        results.append(df)
        if (i // BATCH_SIZE) % 25 == 0:
            elapsed = time.time() - t0
            rate = (i + len(batch)) / max(elapsed, 0.1)
            eta = (n - i - len(batch)) / max(rate, 0.1)
            print(
                f"[uniprot-host] {i + len(batch):,}/{n:,} ({100 * (i + len(batch)) / n:.1f}%) "
                f"rate={rate:.0f}/s eta={eta:.0f}s",
                file=sys.stderr,
            )
        # Polite pause every 50 batches
        if (i // BATCH_SIZE) % 50 == 49:
            time.sleep(1)
    out = pd.concat(results, ignore_index=True)
    out["host_taxid"] = out["host_organism"].map(parse_taxids)
    print(f"[uniprot-host] done: {len(out):,} rows in {time.time() - t0:.0f}s", file=sys.stderr)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_parquet(cache_path, index=False)
        print(f"[uniprot-host] cached → {cache_path}", file=sys.stderr)
    return out


def collect_sample_uniprots(out_dir: Path) -> list[str]:
    """Union of unique uniprots across all out/*/hits_filtered.parquet."""
    uniprots: set[str] = set()
    for p in sorted(out_dir.glob("*/hits_filtered.parquet")):
        df = pd.read_parquet(p, columns=["uniprot_id"])
        uniprots.update(df["uniprot_id"].dropna().unique())
        print(f"[uniprot-host]   {p}: {df['uniprot_id'].nunique():,} uniprots", file=sys.stderr)
    return sorted(uniprots)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument("--out-dir", type=Path, default=Path("/workspace/viral-viz/out"),
                    help="root of sample outputs (used to collect sample-relevant uniprots)")
    ap.add_argument("--refresh", action="store_true",
                    help="Refetch all (including cached) accessions")
    ap.add_argument("--source", choices=["samples", "bfvd"], default="samples",
                    help="samples: union of uniprots in out/*/hits_filtered.parquet (faster); "
                         "bfvd: every uniprot in refs/uniprot_taxid.parquet (heavier, future-proof)")
    args = ap.parse_args()

    cache_path = args.refs_dir / "cache" / "uniprot_host.parquet"

    if args.source == "samples":
        accessions = collect_sample_uniprots(args.out_dir)
    else:
        u = pd.read_parquet(args.refs_dir / "uniprot_taxid.parquet")
        accessions = u["uniprot_id"].dropna().unique().tolist()

    # Filter to valid UniProtKB pattern (drop UPI* and other non-queryable IDs)
    accessions = [a for a in accessions if UNIPROTKB_RE.match(a)]
    print(f"[uniprot-host] candidate UniProtKB accessions: {len(accessions):,}", file=sys.stderr)

    # Incremental: skip ones already in cache (unless --refresh)
    if cache_path.exists() and not args.refresh:
        cached = pd.read_parquet(cache_path)
        already = set(cached["uniprot_id"])
        accessions = [a for a in accessions if a not in already]
        print(f"[uniprot-host]   {len(already):,} cached; fetching {len(accessions):,} new", file=sys.stderr)
    else:
        cached = None

    if not accessions:
        print(f"[uniprot-host] nothing to fetch; cache up-to-date", file=sys.stderr)
        return

    new = fetch_all(accessions, cache_path=None)  # don't overwrite cache yet
    if cached is not None:
        combined = pd.concat([cached, new], ignore_index=True).drop_duplicates(subset=["uniprot_id"], keep="last")
    else:
        combined = new

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(cache_path, index=False)
    n_with_host = combined["host_organism"].astype(bool).sum()
    print(f"[uniprot-host] cache → {cache_path}  ({len(combined):,} total, {n_with_host:,} with host)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
