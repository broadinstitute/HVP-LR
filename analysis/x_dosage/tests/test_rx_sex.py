"""Synthetic-data checks for the (Rx, Ry) karyotype classifier (TODO #5).

These are UNIT TESTS of the classifier math on synthetic mosdepth windows built to
hit target per-copy dosages under the default config. They are NOT real sequencing
data and are NOT a substitute for calibration/validation on labeled controls.

In particular they do NOT validate the long-read aneuploidy classes — those remain
dosage-EXTRAPOLATED (no public long-read aneuploidy data exists; see HANDOFF.md).
Nothing here should be read as evidence that the aneuploidy MEANS are empirically
confirmed on any platform.
"""
from __future__ import annotations

import gzip
import random

import pytest

import rx_sex

CFG = rx_sex.default_config()

# Karyotype (n_X, n_Y) -> expected top call, for the clean euploid + aneuploid grid.
CLEAN_CASES = [
    (2, 0, "46,XX"),
    (1, 1, "46,XY"),
    (1, 0, "45,X"),
    (2, 1, "47,XXY"),
    (1, 2, "47,XYY"),
    (3, 0, "47,XXX"),
    (2, 2, "48,XXYY"),
]


def write_bed(path, rx_target, ry_target, auto_dp=30.0, jitter=0.02, seed=1):
    """Write a synthetic mosdepth regions.bed.gz whose medians give (rx_target, ry_target).

    Autosome windows ~ auto_dp; chrX non-PAR ~ rx_target*auto_dp; chrY euchromatin ~
    ry_target*auto_dp. Mild multiplicative jitter (seeded) keeps it deterministic while
    exercising the bootstrap.
    """
    rng = random.Random(seed)

    def jit(v):
        return v * (1.0 + rng.uniform(-jitter, jitter))

    rows = []
    # 100 autosome windows spread across chr1..chr22
    for i in range(100):
        s = i * 1_000_000
        rows.append((str(1 + i % 22), s, s + 1_000_000, jit(auto_dp)))
    # 60 chrX non-PAR windows (start well past PAR1 end at 2,394,410)
    for i in range(60):
        s = 3_000_000 + i * 1_000_000
        rows.append(("chrX", s, s + 1_000_000, jit(rx_target * auto_dp)))
    # 40 chrY euchromatin windows (inside 2.65-26 Mb)
    for i in range(40):
        s = 3_000_000 + i * 500_000
        rows.append(("chrY", s, s + 500_000, jit(ry_target * auto_dp)))

    with gzip.open(path, "wt") as fh:
        for c, s, e, d in rows:
            fh.write(f"{c}\t{s}\t{e}\t{d}\n")
    return str(path)


def call_bed(path, cfg=CFG):
    """Run the full classify path on a bed; return (ranked_posterior, rx, ry)."""
    random.seed(cfg.seed)  # deterministic bootstrap
    auto, x, y = rx_sex.load(path, cfg)
    _, rx, ry = rx_sex.rx_ry(auto, x, y)
    (_, se_rx, _, se_ry) = rx_sex.boot(x, y, auto, cfg)
    return rx_sex.classify(rx, ry, se_rx, se_ry, cfg), rx, ry


@pytest.mark.parametrize("n_x,n_y,expected", CLEAN_CASES)
def test_clean_karyotypes(tmp_path, n_x, n_y, expected):
    bed = write_bed(
        tmp_path / f"{expected}.bed.gz",
        rx_target=n_x * CFG.rx_unit,
        ry_target=n_y * CFG.ry_unit,
        seed=hash((n_x, n_y)) % 10_000,
    )
    ranked, rx, ry = call_bed(bed)
    top_label, top_p = ranked[0]
    assert top_label == expected, f"{expected}: got {top_label} rx={rx:.3f} ry={ry:.3f}"
    assert top_p > 0.80, f"{expected}: low confidence {top_p:.3f}"


def test_mosaic_goes_to_other(tmp_path):
    # Fractional dosage on both axes (1.5 X copies, 0.5 Y copies) — off every grid point.
    bed = write_bed(tmp_path / "mosaic.bed.gz", rx_target=0.75, ry_target=0.25, seed=7)
    ranked, rx, ry = call_bed(bed)
    top_label, _ = ranked[0]
    assert top_label == "OTHER", f"mosaic should be OTHER, got {top_label} rx={rx:.3f} ry={ry:.3f}"


def test_loy_intermediate_is_ambiguous(tmp_path):
    # Partial loss-of-Y at the XY<->Turner crossover (ry ~0.20, between 46,XY=0.5 and 45,X=0).
    # HANDOFF: LOY males read as an XY<->Turner intermediate — must NOT be a confident call;
    # flag for review. (The ambiguous band is narrow — see the full-LOY test below.)
    bed = write_bed(tmp_path / "loy.bed.gz", rx_target=0.50, ry_target=0.20, seed=11)
    ranked, _, _ = call_bed(bed)
    top_label, top_p = ranked[0]
    assert top_p < 0.90, f"LOY crossover should be ambiguous, got confident {top_label}={top_p:.3f}"
    assert top_label in {"46,XY", "45,X", "OTHER"}, f"unexpected LOY top: {top_label}"


def test_full_loy_reads_as_turner_known_limitation(tmp_path):
    # Documents a KNOWN LIMITATION (HANDOFF): a full-LOY male (ry ~= 0) is confidently
    # absorbed into 45,X — a real male mis-called Turner. The 45,X prior does not prevent
    # this; it is why LOY-range calls must be flagged for human review, not trusted.
    bed = write_bed(tmp_path / "full_loy.bed.gz", rx_target=0.50, ry_target=0.02, seed=13)
    ranked, _, _ = call_bed(bed)
    top_label, _ = ranked[0]
    assert top_label == "45,X"  # not 46,XY — the documented false-Turner failure mode


def test_config_roundtrip_matches_default(tmp_path):
    # config.example.toml is the wired default's serialization: same calls as built-in.
    bed = write_bed(tmp_path / "xy.bed.gz", rx_target=0.5, ry_target=0.5, seed=3)
    file_cfg = rx_sex.load_config("config.example.toml")
    r_default, _, _ = call_bed(bed, CFG)
    r_file, _, _ = call_bed(bed, file_cfg)
    assert r_default[0][0] == r_file[0][0] == "46,XY"
