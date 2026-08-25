"""Unit tests for calibrate.py on SYNTHETIC labeled numbers.

These check the calibration arithmetic (per-copy centers, ref-Y ceiling exclusion,
dosage-linearity ratios, TOML emission) using hand-made (Rx, Ry) values — NOT real
control sequencing. They do not calibrate or validate anything; real calibration
requires the labeled controls in calibration_samples.md. The long-read aneuploidy
classes remain dosage-EXTRAPOLATED regardless of what these tests assert.
"""
from __future__ import annotations

import textwrap

import pytest

import calibrate
import rx_sex

BASE = rx_sex.default_config()
BASE_MAP = rx_sex.DEFAULT_CONFIG


def write_manifest(tmp_path, rows, header="sample_id\tn_x\tn_y\tplatform\try_role\trx\try"):
    p = tmp_path / "controls.tsv"
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return str(p)


def test_units_and_refY_ceiling_excluded(tmp_path):
    rows = [
        "NA12878\t2\t0\tshort_read\tanchor\t1.00\t0.0",
        "HG004\t2\t0\tshort_read\tanchor\t0.98\t0.0",
        "HG003\t1\t1\tshort_read\tanchor\t0.50\t0.45",
        "HG005\t1\t1\tshort_read\tanchor\t0.50\t0.43",
        "HG002\t1\t1\tshort_read\tceiling\t0.50\t0.50",  # ref-Y, must NOT anchor RY_UNIT
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    cal = calibrate.calibrate_platform(controls, "short_read", BASE)
    assert cal.rx_unit == pytest.approx(0.498, abs=1e-3)    # mean of 0.5,0.49,0.5,0.5,0.5
    assert cal.ry_unit == pytest.approx(0.44, abs=1e-3)     # (0.45+0.43)/2, HG002 excluded
    assert cal.ry_ceiling == pytest.approx(0.50, abs=1e-6)  # reported as upper bound only
    assert cal.n_ry_anchor == 2


def test_linearity_pass(tmp_path):
    rows = [
        "XX1\t2\t0\tshort_read\tanchor\t1.00\t0.0",
        "XX2\t2\t0\tshort_read\tanchor\t0.99\t0.0",
        "XY1\t1\t1\tshort_read\tanchor\t0.50\t0.45",
        "XY2\t1\t1\tshort_read\tanchor\t0.49\t0.44",
        "XYY1\t1\t2\tshort_read\tanchor\t0.50\t0.89",  # 2Y ~ 2*1Y
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    rungs = {r.axis.split()[0]: r for r in calibrate.linearity_check(controls, "short_read")}
    assert rungs["X"].status == "PASS"
    assert rungs["Y"].status == "PASS"
    assert rungs["Y"].ratio == pytest.approx(0.89 / (2 * 0.445), abs=1e-3)


def test_linearity_insufficient_without_2Y(tmp_path):
    rows = [
        "XX1\t2\t0\tshort_read\tanchor\t1.00\t0.0",
        "XY1\t1\t1\tshort_read\tanchor\t0.50\t0.45",
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    rungs = {r.axis.split()[0]: r for r in calibrate.linearity_check(controls, "short_read")}
    assert rungs["X"].status == "PASS"          # XX vs 2*XY present
    assert rungs["Y"].status == "INSUFFICIENT"  # no XYY -> not fabricated


def test_fallback_when_only_ceiling_male(tmp_path):
    # Only a ref-Y ceiling male, no anchors: RY_UNIT must fall back to base, with a warning
    # (never anchor RY_UNIT on the mappability ceiling).
    rows = [
        "XX1\t2\t0\tshort_read\tanchor\t1.00\t0.0",
        "HG002\t1\t1\tshort_read\tceiling\t0.50\t0.50",
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    cal = calibrate.calibrate_platform(controls, "short_read", BASE)
    assert cal.ry_unit == pytest.approx(BASE.ry_unit)
    assert cal.n_ry_anchor == 0
    assert any("anchor" in w.lower() for w in cal.warnings)


def test_long_read_extrapolation_caveat(tmp_path):
    rows = [
        "HG001\t2\t0\tlong_read\tanchor\t1.00\t0.0",
        "HG005\t1\t1\tlong_read\tanchor\t0.50\t0.40",
        "HG006\t1\t1\tlong_read\tanchor\t0.49\t0.39",
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    cal = calibrate.calibrate_platform(controls, "long_read", BASE)
    assert any("EXTRAPOLATED" in w for w in cal.warnings), "LR aneuploidy caveat must persist"


def test_emit_toml_roundtrips(tmp_path):
    rows = [
        "XX1\t2\t0\tshort_read\tanchor\t1.00\t0.0",
        "XY1\t1\t1\tshort_read\tanchor\t0.50\t0.45",
        "XY2\t1\t1\tshort_read\tanchor\t0.50\t0.43",
    ]
    controls = calibrate.read_manifest(write_manifest(tmp_path, rows), BASE)
    cal = calibrate.calibrate_platform(controls, "short_read", BASE)
    toml_text = calibrate.emit_toml(cal, BASE_MAP, "controls.tsv")
    out = tmp_path / "short_read.toml"
    out.write_text(toml_text + "\n")

    cfg = rx_sex.load_config(str(out))          # must parse as valid TOML
    assert cfg.platform == "short_read"
    assert cfg.rx_unit == pytest.approx(cal.rx_unit, abs=1e-4)
    assert cfg.ry_unit == pytest.approx(cal.ry_unit, abs=1e-4)
    assert len(cfg.karyotypes) == len(BASE.karyotypes)


def test_manifest_requires_signal_columns(tmp_path):
    # Missing both 'regions' and rx/ry -> hard error, not a silent empty calibration.
    bad = tmp_path / "bad.tsv"
    bad.write_text(textwrap.dedent("""\
        sample_id\tn_x\tn_y\tplatform
        XX1\t2\t0\tshort_read
    """))
    with pytest.raises(ValueError, match="regions.*rx.*ry|rx.*ry"):
        calibrate.read_manifest(str(bad), BASE)
