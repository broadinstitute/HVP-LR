#!/usr/bin/env python3
"""Emit a per-platform rx_sex config from labeled controls (TODO #3).

Reads a manifest of known-karyotype control samples, measures each sample's
(Rx, Ry) with the same math rx_sex.py uses, and derives the calibration constants:

  RX_UNIT        = per-copy X center  = mean over samples of Rx / n_X   (n_X > 0)
  RY_UNIT        = per-copy Y center  = mean over ANCHOR males of Ry / n_Y (n_Y > 0)
  SIGMA_RX_MODEL = within-group spread = SD of residuals Rx - n_X*RX_UNIT
  SIGMA_RY_MODEL = within-group spread = SD of residuals Ry - n_Y*RY_UNIT  (anchors)

One config is emitted per `platform` value in the manifest (short-read and
long-read have different Y mappability, so never share constants).

RY_UNIT de-bias: HG002 maps to its OWN reference Y (mappability ceiling), so it
biases RY_UNIT high. Mark it (and any ref-Y sample) `ry_role = ceiling` in the
manifest; it is excluded from the RY_UNIT/SIGMA_RY fit and reported as an upper
bound only. Anchor RY_UNIT on non-HG002 males (HG005/HG006, 1kGP males).

Dosage-linearity check (does the 2-fold hold?):
  X:  2X_center vs 2 * 1X_center     (e.g. XX vs 2*XY)
  Y:  2Y_center vs 2 * 1Y_center      (e.g. XYY vs 2*XY)
Reported as a ratio with PASS/FAIL; if a rung has no controls it is INSUFFICIENT
(never fabricated).

CAVEAT (do not remove): only euploid per-copy UNITS are directly calibratable on
the long-read path (GIAB euploids). Long-read ANEUPLOIDY class means stay
dosage-EXTRAPOLATED — there is no public long-read aneuploidy data to fit or
validate them. This tool does not manufacture that data.

Manifest: TSV with a header. Required columns: sample_id, n_x, n_y, platform.
Per-sample signal comes from EITHER `regions` (a mosdepth *.regions.bed(.gz) path,
measured here) OR precomputed `rx` and `ry` columns. Optional: ry_role
(anchor|ceiling, default anchor).

Stdlib only (tomllib to read the base config; a small hand-written TOML emitter).
Usage:
  python calibrate.py --manifest controls.tsv [--base-config config.example.toml]
                      [--out-dir configs] [--linearity-tol 0.15]
"""
from __future__ import annotations

import argparse
import copy
import csv
import statistics as st
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

import rx_sex

LINEARITY_TOL = 0.15  # |ratio - 1| tolerance for the 2-fold dosage check


@dataclass
class Control:
    sample_id: str
    n_x: int
    n_y: int
    platform: str
    ry_role: str  # "anchor" | "ceiling"
    rx: float
    ry: float


def _measure(regions: str, cfg: rx_sex.Config) -> tuple[float, float]:
    """Point (Rx, Ry) for one mosdepth regions file (no bootstrap needed for calibration)."""
    auto, x, y = rx_sex.load(regions, cfg)
    if not auto or not x:
        raise ValueError(f"{regions}: insufficient autosome/chrX windows")
    _, rx, ry = rx_sex.rx_ry(auto, x, y)
    return rx, ry


def read_manifest(path: str, base_cfg: rx_sex.Config) -> list[Control]:
    controls: list[Control] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = set(reader.fieldnames or [])
        for req in ("sample_id", "n_x", "n_y", "platform"):
            if req not in cols:
                raise ValueError(f"manifest missing required column: {req}")
        has_precomputed = {"rx", "ry"} <= cols
        has_regions = "regions" in cols
        if not (has_precomputed or has_regions):
            raise ValueError("manifest needs either 'regions' or both 'rx' and 'ry' columns")
        for row in reader:
            sid = row["sample_id"].strip()
            if not sid or sid.startswith("#"):
                continue
            n_x, n_y = int(row["n_x"]), int(row["n_y"])
            regions = (row.get("regions") or "").strip()
            if regions:
                rx, ry = _measure(regions, base_cfg)
            elif row.get("rx", "") != "" and row.get("ry", "") != "":
                rx, ry = float(row["rx"]), float(row["ry"])
            else:
                raise ValueError(f"{sid}: no regions path and no rx/ry values")
            controls.append(
                Control(
                    sample_id=sid,
                    n_x=n_x,
                    n_y=n_y,
                    platform=row["platform"].strip(),
                    ry_role=(row.get("ry_role") or "anchor").strip().lower(),
                    rx=rx,
                    ry=ry,
                )
            )
    if not controls:
        raise ValueError("manifest has no usable rows")
    return controls


@dataclass
class Calibration:
    platform: str
    rx_unit: float
    ry_unit: float
    sigma_rx_model: float
    sigma_ry_model: float
    ry_ceiling: float | None          # HG002-style upper bound, informational
    n_rx: int
    n_ry_anchor: int
    warnings: list[str]


def calibrate_platform(
    controls: list[Control], platform: str, base_cfg: rx_sex.Config
) -> Calibration:
    warns: list[str] = []
    rows = [c for c in controls if c.platform == platform]

    # RX_UNIT from every sample carrying X (all karyotypes have n_x >= 1 here).
    x_rows = [c for c in rows if c.n_x > 0]
    if not x_rows:
        raise ValueError(f"[{platform}] no samples with n_x>0 to fit RX_UNIT")
    rx_unit = st.fmean([c.rx / c.n_x for c in x_rows])

    # RY_UNIT from ANCHOR males only; ref-Y ceilings excluded (reported separately).
    y_anchor = [c for c in rows if c.n_y > 0 and c.ry_role == "anchor"]
    y_ceiling = [c for c in rows if c.n_y > 0 and c.ry_role == "ceiling"]
    ry_ceiling_val: float | None = None
    if y_ceiling:
        ry_ceiling_val = max(c.ry / c.n_y for c in y_ceiling)
    if y_anchor:
        ry_unit = st.fmean([c.ry / c.n_y for c in y_anchor])
    else:
        ry_unit = base_cfg.ry_unit
        warns.append(
            "no ANCHOR males (ry_role=anchor, n_y>0) — RY_UNIT kept at base default; "
            "do NOT anchor on a ref-Y ceiling sample"
        )

    # Dispersions = SD of residuals about the fitted per-copy line.
    if len(x_rows) >= 2:
        sigma_rx = st.pstdev([c.rx - c.n_x * rx_unit for c in x_rows])
    else:
        sigma_rx = base_cfg.sigma_rx_model
        warns.append("<2 X samples — SIGMA_RX_MODEL kept at base default")
    if len(y_anchor) >= 2:
        sigma_ry = st.pstdev([c.ry - c.n_y * ry_unit for c in y_anchor])
    else:
        sigma_ry = base_cfg.sigma_ry_model
        warns.append("<2 anchor males — SIGMA_RY_MODEL kept at base default")

    if platform == "long_read":
        warns.append(
            "long_read: only euploid UNITS calibrated; aneuploidy class means remain "
            "dosage-EXTRAPOLATED (no public long-read aneuploidy data)"
        )

    return Calibration(
        platform=platform,
        rx_unit=rx_unit,
        ry_unit=ry_unit,
        sigma_rx_model=sigma_rx,
        sigma_ry_model=sigma_ry,
        ry_ceiling=ry_ceiling_val,
        n_rx=len(x_rows),
        n_ry_anchor=len(y_anchor),
        warnings=warns,
    )


@dataclass
class Rung:
    axis: str
    one_center: float | None
    two_center: float | None
    ratio: float | None
    status: str  # PASS | FAIL | INSUFFICIENT


def _center(vals: list[float]) -> float | None:
    return st.fmean(vals) if vals else None


def linearity_check(
    controls: list[Control], platform: str, tol: float = LINEARITY_TOL
) -> list[Rung]:
    """2X ≈ 2·1X on the X axis; 2Y ≈ 2·1Y on the Y axis (anchors only for Y)."""
    rows = [c for c in controls if c.platform == platform]
    out: list[Rung] = []

    one_x = _center([c.rx for c in rows if c.n_x == 1])
    two_x = _center([c.rx for c in rows if c.n_x == 2])
    out.append(_rung("X (2X vs 2*1X)", one_x, two_x, tol))

    anchors = [c for c in rows if c.ry_role == "anchor"]
    one_y = _center([c.ry for c in anchors if c.n_y == 1])
    two_y = _center([c.ry for c in anchors if c.n_y == 2])
    out.append(_rung("Y (2Y vs 2*1Y)", one_y, two_y, tol))
    return out


def _rung(axis: str, one: float | None, two: float | None, tol: float) -> Rung:
    if one is None or two is None or one == 0:
        return Rung(axis, one, two, None, "INSUFFICIENT")
    ratio = two / (2.0 * one)
    return Rung(axis, one, two, ratio, "PASS" if abs(ratio - 1.0) <= tol else "FAIL")


def _fmt_karyos(karyos: list[dict]) -> str:
    blocks = []
    for k in karyos:
        label = k["label"]
        comment = "    # Klinefelter" if label == "47,XXY" else (
            "     # Turner" if label == "45,X" else "")
        blocks.append(
            f'[[karyotypes]]\nlabel = "{label}"{comment}\n'
            f'n_x = {int(k["n_x"])}\nn_y = {int(k["n_y"])}\nprior = {float(k["prior"])}\n'
        )
    return "\n".join(blocks)


def emit_toml(cal: Calibration, base_map: dict, manifest: str) -> str:
    """Render a valid TOML config: base regions/other/bootstrap/karyotypes, calibrated units."""
    reg = base_map["regions"]
    other = base_map["other"]
    boot_cfg = base_map.get("bootstrap", {"n_boot": 2000, "seed": 0})
    ceiling = (
        f"  # ANCHOR fit; ref-Y ceiling (excluded) = {cal.ry_ceiling:.4f}"
        if cal.ry_ceiling is not None
        else ""
    )
    header = [
        f"# Calibrated by calibrate.py from manifest: {manifest}",
        f"# platform={cal.platform}  n_X_samples={cal.n_rx}  n_Y_anchor={cal.n_ry_anchor}",
    ]
    if cal.platform == "long_read":
        header.append(
            "# NOTE: euploid units only; aneuploidy class means are dosage-EXTRAPOLATED "
            "(no public long-read aneuploidy data)."
        )
    for w in cal.warnings:
        header.append(f"# WARN: {w}")
    return (
        "\n".join(header)
        + "\n\n"
        + f'platform = "{cal.platform}"\n'
        + f'reference = "{base_map.get("reference", "chm13v2.0_maskedY.rCRS")}"\n\n'
        + "[regions]\n"
        + f"par1          = [{reg['par1'][0]}, {reg['par1'][1]}]\n"
        + f"par2          = [{reg['par2'][0]}, {reg['par2'][1]}]\n"
        + f"y_euchromatin = [{reg['y_euchromatin'][0]}, {reg['y_euchromatin'][1]}]\n\n"
        + "[units]\n"
        + f"rx_unit = {cal.rx_unit:.4f}\n"
        + f"ry_unit = {cal.ry_unit:.4f}{ceiling}\n\n"
        + "[dispersion]\n"
        + f"sigma_rx_model = {cal.sigma_rx_model:.4f}\n"
        + f"sigma_ry_model = {cal.sigma_ry_model:.4f}\n\n"
        + "[other]\n"
        + f"rx_max      = {float(other['rx_max'])}\n"
        + f"ry_max      = {float(other['ry_max'])}\n"
        + f"prior_other = {float(other['prior_other'])}\n\n"
        + "[bootstrap]\n"
        + f"n_boot = {int(boot_cfg.get('n_boot', 2000))}\n"
        + f"seed   = {int(boot_cfg.get('seed', 0))}\n\n"
        + _fmt_karyos(base_map["karyotypes"])
    )


def load_base_map(path: str | None) -> dict:
    if path:
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    return copy.deepcopy(rx_sex.DEFAULT_CONFIG)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", required=True, help="TSV of labeled controls")
    ap.add_argument("--base-config", help="base TOML (regions/other/karyotypes); default built-in")
    ap.add_argument("--out-dir", default="configs", help="where to write <platform>.toml")
    ap.add_argument("--linearity-tol", type=float, default=LINEARITY_TOL)
    args = ap.parse_args(argv)

    base_map = load_base_map(args.base_config)
    base_cfg = rx_sex.config_from_mapping(base_map)
    controls = read_manifest(args.manifest, base_cfg)

    platforms = sorted({c.platform for c in controls})
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for platform in platforms:
        cal = calibrate_platform(controls, platform, base_cfg)
        toml_text = emit_toml(cal, base_map, args.manifest)
        out_path = out_dir / f"{platform}.toml"
        out_path.write_text(toml_text + "\n")

        print(f"\n===== platform: {platform} =====")
        print(
            f"RX_UNIT={cal.rx_unit:.4f}  RY_UNIT={cal.ry_unit:.4f}  "
            f"SIGMA_RX={cal.sigma_rx_model:.4f}  SIGMA_RY={cal.sigma_ry_model:.4f}"
        )
        if cal.ry_ceiling is not None:
            print(f"  (ref-Y ceiling, excluded from RY_UNIT: {cal.ry_ceiling:.4f})")
        print(f"  fit from n_X={cal.n_rx} X-bearing, n_Y_anchor={cal.n_ry_anchor} anchor males")
        for w in cal.warnings:
            print(f"  WARN: {w}")

        print("  dosage linearity:")
        for r in linearity_check(controls, platform, args.linearity_tol):
            if r.status == "INSUFFICIENT":
                print(f"    {r.axis:20s} INSUFFICIENT (no controls on one rung)")
            else:
                print(
                    f"    {r.axis:20s} 1x={r.one_center:.4f} 2x={r.two_center:.4f} "
                    f"ratio={r.ratio:.3f} -> {r.status}"
                )
        print(f"  wrote {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
