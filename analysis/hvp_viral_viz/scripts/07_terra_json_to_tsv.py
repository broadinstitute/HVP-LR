"""Convert Terra get_entities JSON dumps to TSV.

Input shape: {"entity_type": str, "count": int, "entities": [
    {"name": str, "entityType": str, "attributes": {col: value, ...}}, ...
]}

Attribute values may be:
- scalar (str/int/float/bool/null) → written as-is
- {"itemsType": "AttributeValue"|"EntityReference", "items": [...]} → joined
  with ";" (entityName extracted for EntityReference)
- dict otherwise → JSON-encoded

One TSV per input JSON. Column order: name, then attribute keys sorted.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def flatten_value(v):
    if v is None:
        return ""
    if isinstance(v, (str, int, float, bool)):
        return str(v)
    if isinstance(v, dict):
        if "itemsType" in v and "items" in v:
            items = v["items"]
            if v["itemsType"] == "EntityReference":
                return ";".join(str(x.get("entityName", "")) for x in items)
            return ";".join(str(x) for x in items)
        return json.dumps(v, separators=(",", ":"))
    if isinstance(v, list):
        return ";".join(str(x) for x in v)
    return str(v)


def convert(json_path: Path, tsv_path: Path) -> tuple[int, int]:
    payload = json.loads(json_path.read_text())
    entities = payload["entities"]

    # Union of attribute keys across all rows (deterministic order).
    keys = sorted({k for e in entities for k in e.get("attributes", {}).keys()})
    cols = ["name"] + keys

    with tsv_path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in entities:
            attrs = e.get("attributes", {})
            row = [e.get("name", "")] + [flatten_value(attrs.get(k)) for k in keys]
            # Strip tabs/newlines inside cells to keep TSV well-formed.
            row = [c.replace("\t", " ").replace("\n", " ").replace("\r", " ") for c in row]
            fh.write("\t".join(row) + "\n")
    return len(entities), len(cols)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, type=Path,
                    help="Directory containing <table>.json files to convert in place to <table>.tsv")
    args = ap.parse_args()

    json_paths = sorted(args.dir.glob("*.json"))
    if not json_paths:
        print(f"[convert] no *.json in {args.dir}", file=sys.stderr)
        sys.exit(1)

    for jp in json_paths:
        tsv = jp.with_suffix(".tsv")
        n_rows, n_cols = convert(jp, tsv)
        print(f"[convert] {jp.name} → {tsv.name}  ({n_rows} rows × {n_cols} cols)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
