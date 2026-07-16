"""Mirror of the inline Python heredoc in
wdl/tasks/ProteinAnnotation/ProteinAnnotationHelpers.wdl :: AnnotationTransfer.

Functions are factored out for unit testing. If you change one, change the
heredoc to match (or vice versa) and re-run miniwdl --strict + pytest.
"""
from __future__ import annotations

import csv
import os
from typing import Iterable


class SchemaError(ValueError):
    pass


def validate_hits_columns(hits_cols: list[str]) -> tuple[int, int | None]:
    """Return (evalue_idx, bits_idx_or_None). Raise SchemaError on bad schema."""
    if len(hits_cols) < 3 or hits_cols[0] != "query" or hits_cols[1] != "target":
        raise SchemaError(
            f"hits_columns must start with ['query', 'target', ...]; got {hits_cols}"
        )
    if "evalue" not in hits_cols:
        raise SchemaError(f"hits_columns must include 'evalue'; got {hits_cols}")
    evalue_idx = hits_cols.index("evalue")
    bits_idx = hits_cols.index("bits") if "bits" in hits_cols else None
    return evalue_idx, bits_idx


def reduce_best_hits(
    hits_path: str | os.PathLike,
    hits_cols: list[str],
) -> dict[str, list[str]]:
    """One row per query — best foldseek hit by lowest e-value, tiebreaker higher bits.

    Rows shorter than len(hits_cols) and rows with non-float evalue are skipped.
    Returns {query -> row[:ncol]}.
    """
    evalue_idx, bits_idx = validate_hits_columns(hits_cols)
    ncol = len(hits_cols)
    best: dict[str, tuple[tuple[float, float], list[str]]] = {}
    with open(hits_path, newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or len(row) < ncol:
                continue
            q = row[0]
            try:
                ev_f = float(row[evalue_idx])
            except ValueError:
                continue
            if bits_idx is not None:
                try:
                    bt_f = float(row[bits_idx])
                except ValueError:
                    bt_f = 0.0
                key = (ev_f, -bt_f)
            else:
                key = (ev_f, 0.0)
            if q not in best or key < best[q][0]:
                best[q] = (key, row[:ncol])
    return {q: v[1] for q, v in best.items()}


def load_metadata(
    meta_path: str | os.PathLike,
    has_header: bool,
    explicit_cols: list[str] | None,
) -> tuple[list[str] | None, dict[str, list[str]], int, int]:
    """Load reference metadata TSV.

    Returns (header_or_None, key->row, n_meta_rows, n_dup_keys).
    Duplicate keys collapse to last-row-wins.
    Synthetic 'meta_col_N' headers are generated only when has_header is False
    AND explicit_cols is falsy; they are seeded from the first data row's width.
    """
    meta: dict[str, list[str]] = {}
    header: list[str] | None = None
    n_rows = 0
    n_dups = 0
    with open(meta_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        if has_header:
            header = next(reader, None)
        for r in reader:
            if not r:
                continue
            n_rows += 1
            if header is None and not explicit_cols:
                header = [f"meta_col_{i+1}" for i in range(len(r))]
            if r[0] in meta:
                n_dups += 1
            meta[r[0]] = r
    if not has_header and explicit_cols:
        header = list(explicit_cols)
    return header, meta, n_rows, n_dups


def join_annotations(
    best: dict[str, list[str]],
    hits_cols: list[str],
    meta_header: list[str] | None,
    meta: dict[str, list[str]],
) -> tuple[list[str], list[list[str]], int]:
    """Left-join best hits against metadata.

    Returns (out_header, rows, n_annotated). When meta_header is None (e.g.
    metadata file was empty), the output is identical-shape to best_hits.
    """
    if meta_header is None:
        out_header = list(hits_cols)
        rows = [list(best[q]) for q in sorted(best)]
        return out_header, rows, 0

    ann_cols = ["target_meta_" + c for c in meta_header]
    out_header = list(hits_cols) + ann_cols
    rows: list[list[str]] = []
    n_annotated = 0
    for q in sorted(best):
        row = list(best[q])
        target = row[1]
        m = meta.get(target)
        if m is not None:
            m_padded = (m + [""] * len(ann_cols))[: len(ann_cols)]
            row += m_padded
            n_annotated += 1
        else:
            row += [""] * len(ann_cols)
        rows.append(row)
    return out_header, rows, n_annotated


def write_tsv(path: str | os.PathLike, header: Iterable[str], rows: Iterable[Iterable[str]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(list(header))
        for r in rows:
            w.writerow(list(r))
