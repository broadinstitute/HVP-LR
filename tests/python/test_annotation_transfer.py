"""pytest for the AnnotationTransfer join logic.

Run from repo root: `python3 -m pytest tests/python -v`

Tests target tests/python/annotation_transfer_core.py, which mirrors the
Python heredoc inside the AnnotationTransfer task in
wdl/tasks/ProteinAnnotation/ProteinAnnotationHelpers.wdl. The two MUST stay
in sync.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from annotation_transfer_core import (
    SchemaError,
    join_annotations,
    load_metadata,
    reduce_best_hits,
    validate_hits_columns,
    write_tsv,
)


DEFAULT_COLS = ["query", "target", "evalue", "bits", "fident", "alnlen", "prob"]


def _write_hits(path: Path, rows: list[list[str]]) -> None:
    path.write_text("".join("\t".join(r) + "\n" for r in rows))


# ---- validate_hits_columns ----

def test_validate_default_ok():
    ev, bt = validate_hits_columns(DEFAULT_COLS)
    assert ev == 2 and bt == 3


def test_validate_no_bits():
    ev, bt = validate_hits_columns(["query", "target", "evalue"])
    assert ev == 2 and bt is None


def test_validate_missing_query_first():
    with pytest.raises(SchemaError):
        validate_hits_columns(["q", "target", "evalue"])


def test_validate_missing_evalue():
    with pytest.raises(SchemaError):
        validate_hits_columns(["query", "target", "bits"])


def test_validate_too_short():
    with pytest.raises(SchemaError):
        validate_hits_columns(["query", "target"])


# ---- reduce_best_hits ----

def test_reduce_empty(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text("")
    assert reduce_best_hits(p, DEFAULT_COLS) == {}


def test_reduce_single_per_query(tmp_path):
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t1", "1e-10", "200", "0.9", "100", "0.95"],
        ["q2", "t2", "1e-5",  "150", "0.8", "120", "0.90"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    assert set(best.keys()) == {"q1", "q2"}
    assert best["q1"][1] == "t1"


def test_reduce_multi_lowest_evalue_wins(tmp_path):
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t_lo",  "1e-2", "100", "x", "x", "x"],
        ["q1", "t_med", "1e-8", "200", "x", "x", "x"],
        ["q1", "t_hi",  "1e-3", "150", "x", "x", "x"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    assert best["q1"][1] == "t_med"


def test_reduce_evalue_tie_higher_bits_wins(tmp_path):
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t_low_bits",  "1e-10", "100", "x", "x", "x"],
        ["q1", "t_high_bits", "1e-10", "300", "x", "x", "x"],
        ["q1", "t_mid_bits",  "1e-10", "200", "x", "x", "x"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    assert best["q1"][1] == "t_high_bits"


def test_reduce_evalue_tie_no_bits_column_first_wins(tmp_path):
    cols = ["query", "target", "evalue"]
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t_first",  "1e-10"],
        ["q1", "t_second", "1e-10"],
    ])
    best = reduce_best_hits(p, cols)
    assert best["q1"][1] == "t_first"  # first-wins on strict-less-than insert


def test_reduce_skips_short_rows(tmp_path):
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t1"],
        ["q2", "t2", "1e-5", "150", "x", "x", "x"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    assert set(best.keys()) == {"q2"}


def test_reduce_skips_nonfloat_evalue(tmp_path):
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t1", "bogus", "100", "x", "x", "x"],
        ["q1", "t2", "1e-5",  "200", "x", "x", "x"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    assert best["q1"][1] == "t2"


def test_reduce_nonfloat_bits_falls_back_to_zero(tmp_path):
    """bits is the tie-break; non-float bits should not crash, treated as 0."""
    p = tmp_path / "hits.tsv"
    _write_hits(p, [
        ["q1", "t_clean", "1e-10", "100",  "x", "x", "x"],
        ["q1", "t_dirty", "1e-10", "junk", "x", "x", "x"],
    ])
    best = reduce_best_hits(p, DEFAULT_COLS)
    # t_clean bits=100 > 0; t_clean wins
    assert best["q1"][1] == "t_clean"


# ---- load_metadata ----

def test_load_metadata_with_header(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("acc\torg\tfunc\n"
                 "U1\tBacillus\tlysozyme\n"
                 "U2\tE.coli\tholin\n")
    header, meta, n_rows, n_dups = load_metadata(p, has_header=True, explicit_cols=None)
    assert header == ["acc", "org", "func"]
    assert meta["U1"][1] == "Bacillus"
    assert n_rows == 2 and n_dups == 0


def test_load_metadata_duplicate_keys_last_wins(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("acc\torg\n"
                 "U1\tFirst\n"
                 "U1\tSecond\n"
                 "U1\tThird\n")
    header, meta, n_rows, n_dups = load_metadata(p, has_header=True, explicit_cols=None)
    assert meta["U1"][1] == "Third"
    assert n_rows == 3 and n_dups == 2


def test_load_metadata_no_header_synthetic(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("U1\tFoo\tBar\n"
                 "U2\tBaz\tQux\n")
    header, meta, n_rows, n_dups = load_metadata(p, has_header=False, explicit_cols=None)
    assert header == ["meta_col_1", "meta_col_2", "meta_col_3"]
    assert meta["U1"] == ["U1", "Foo", "Bar"]


def test_load_metadata_no_header_explicit_cols(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("U1\tFoo\tBar\n")
    header, _, _, _ = load_metadata(p, has_header=False, explicit_cols=["acc", "a", "b"])
    assert header == ["acc", "a", "b"]


def test_load_metadata_empty_file(tmp_path):
    p = tmp_path / "meta.tsv"
    p.write_text("")
    header, meta, n_rows, n_dups = load_metadata(p, has_header=True, explicit_cols=None)
    assert header is None
    assert meta == {}
    assert n_rows == 0 and n_dups == 0


# ---- join_annotations ----

def test_join_no_metadata_passthrough():
    best = {"q1": ["q1", "t1", "1e-10", "100", "x", "x", "x"]}
    out_h, rows, n_ann = join_annotations(best, DEFAULT_COLS, meta_header=None, meta={})
    assert out_h == DEFAULT_COLS
    assert rows == [["q1", "t1", "1e-10", "100", "x", "x", "x"]]
    assert n_ann == 0


def test_join_hit_target_in_meta():
    best = {"q1": ["q1", "U1", "1e-10", "100", "x", "x", "x"]}
    out_h, rows, n_ann = join_annotations(
        best,
        DEFAULT_COLS,
        meta_header=["acc", "org", "func"],
        meta={"U1": ["U1", "Bacillus", "lysozyme"]},
    )
    assert out_h == DEFAULT_COLS + ["target_meta_acc", "target_meta_org", "target_meta_func"]
    assert rows[0][-3:] == ["U1", "Bacillus", "lysozyme"]
    assert n_ann == 1


def test_join_hit_target_not_in_meta_empty_fill():
    best = {"q1": ["q1", "Uxxx", "1e-10", "100", "x", "x", "x"]}
    out_h, rows, n_ann = join_annotations(
        best,
        DEFAULT_COLS,
        meta_header=["acc", "org", "func"],
        meta={"U1": ["U1", "Bacillus", "lysozyme"]},
    )
    assert rows[0][-3:] == ["", "", ""]
    assert n_ann == 0


def test_join_meta_row_short_padded():
    """Metadata row with fewer columns than header → trailing empties."""
    best = {"q1": ["q1", "U1", "1e-10", "100", "x", "x", "x"]}
    out_h, rows, n_ann = join_annotations(
        best,
        DEFAULT_COLS,
        meta_header=["acc", "org", "func"],
        meta={"U1": ["U1", "Bacillus"]},  # missing 'func'
    )
    assert rows[0][-3:] == ["U1", "Bacillus", ""]
    assert n_ann == 1


def test_join_output_sorted_by_query():
    best = {
        "q_z": ["q_z", "tZ", "1e-1", "10", "x", "x", "x"],
        "q_a": ["q_a", "tA", "1e-2", "20", "x", "x", "x"],
        "q_m": ["q_m", "tM", "1e-3", "30", "x", "x", "x"],
    }
    _, rows, _ = join_annotations(best, DEFAULT_COLS, meta_header=None, meta={})
    assert [r[0] for r in rows] == ["q_a", "q_m", "q_z"]


# ---- write_tsv smoke ----

def test_write_tsv_roundtrip(tmp_path):
    p = tmp_path / "out.tsv"
    write_tsv(p, ["a", "b"], [["1", "2"], ["3", "4"]])
    assert p.read_text() == "a\tb\n1\t2\n3\t4\n"
