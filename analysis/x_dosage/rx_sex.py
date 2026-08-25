#!/usr/bin/env python3
"""Sex-chromosome karyotyping from mosdepth windows via joint (Rx, Ry) posterior.

Input: mosdepth `*.regions.bed.gz` (run with --by 1000000 -Q 20) on T2T-CHM13v2.0
(chm13v2.0_maskedY.rCRS.fa recommended).

  Rx = median(chrX non-PAR windows)      / median(autosome windows)   ~ n_X * RX_UNIT
  Ry = median(chrY euchromatin windows)  / median(autosome windows)   ~ n_Y * RY_UNIT

Each karyotype h = (n_X, n_Y) predicts (mu_rx, mu_ry). Confidence is the posterior
over the karyotype set + a diffuse OTHER class (mosaics/contamination/structural).
Per-sample sigma = sqrt(SIGMA_MODEL^2 + bootstrap_SE^2), independently for Rx and Ry.

CALIBRATE on known controls before trusting magnitudes — especially RY_UNIT and
SIGMA_RY_MODEL: Y mappability and loss-of-Y make the real per-copy Ry < 0.5 and its
spread wider. Y-only-distinguished pairs (XX/XXY, XY/XYY) live or die on that.

The per-copy units, dispersions, region coordinates and karyotype table now come
from a TOML config (one per platform; format in config.example.toml). With no
--config the built-in DEFAULT_CONFIG (the historical hardcoded values) is used, so
the CLI still runs out of the box. Calibrate real numbers with calibrate.py.

Coordinates: T2T-CHM13v2.0. Core is stdlib only (tomllib is stdlib on 3.11+).
Usage: python rx_sex.py [--config platform.toml] A.regions.bed.gz [B ...]
"""
from __future__ import annotations

import argparse
import gzip
import math
import random
import statistics as st
import tomllib
from dataclasses import dataclass
from pathlib import Path

AUTO = {str(i) for i in range(1, 23)}

# Built-in default = the historical hardcoded constants, in the config.example.toml
# shape so the file path and the default path share one builder. CALIBRATE before
# trusting magnitudes (see calibrate.py).
DEFAULT_CONFIG: dict = {
    "platform": "default",
    "reference": "chm13v2.0_maskedY.rCRS",
    "regions": {
        "par1": [0, 2_394_410],
        "par2": [153_925_834, 154_259_566],
        "y_euchromatin": [2_650_000, 26_000_000],
    },
    "units": {"rx_unit": 0.50, "ry_unit": 0.50},
    "dispersion": {"sigma_rx_model": 0.05, "sigma_ry_model": 0.08},
    "other": {"rx_max": 2.5, "ry_max": 1.5, "prior_other": 0.05},
    "bootstrap": {"n_boot": 2000, "seed": 0},
    "karyotypes": [
        {"label": "46,XX", "n_x": 2, "n_y": 0, "prior": 1.00},
        {"label": "46,XY", "n_x": 1, "n_y": 1, "prior": 1.00},
        {"label": "47,XXY", "n_x": 2, "n_y": 1, "prior": 0.10},   # Klinefelter
        {"label": "47,XYY", "n_x": 1, "n_y": 2, "prior": 0.10},
        {"label": "47,XXX", "n_x": 3, "n_y": 0, "prior": 0.10},
        {"label": "45,X", "n_x": 1, "n_y": 0, "prior": 0.05},     # Turner
        {"label": "48,XXYY", "n_x": 2, "n_y": 2, "prior": 0.01},
    ],
}


@dataclass(frozen=True)
class Karyotype:
    label: str
    n_x: int
    n_y: int
    prior: float


@dataclass(frozen=True)
class Config:
    platform: str
    reference: str
    par1: tuple[int, int]
    par2: tuple[int, int]
    y_euchromatin: tuple[int, int]
    rx_unit: float
    ry_unit: float
    sigma_rx_model: float
    sigma_ry_model: float
    rx_max: float
    ry_max: float
    prior_other: float
    n_boot: int
    seed: int
    karyotypes: tuple[Karyotype, ...]


def _pair(v) -> tuple[int, int]:
    a, b = v
    return (int(a), int(b))


def config_from_mapping(d: dict) -> Config:
    """Build a Config from a parsed-TOML-shaped mapping. Shared by default + file."""
    regions = d["regions"]
    units = d["units"]
    disp = d["dispersion"]
    other = d["other"]
    boot_cfg = d.get("bootstrap", {})
    karyos = tuple(
        Karyotype(k["label"], int(k["n_x"]), int(k["n_y"]), float(k["prior"]))
        for k in d["karyotypes"]
    )
    if not karyos:
        raise ValueError("config has no karyotypes")
    return Config(
        platform=str(d.get("platform", "unknown")),
        reference=str(d.get("reference", "")),
        par1=_pair(regions["par1"]),
        par2=_pair(regions["par2"]),
        y_euchromatin=_pair(regions["y_euchromatin"]),
        rx_unit=float(units["rx_unit"]),
        ry_unit=float(units["ry_unit"]),
        sigma_rx_model=float(disp["sigma_rx_model"]),
        sigma_ry_model=float(disp["sigma_ry_model"]),
        rx_max=float(other["rx_max"]),
        ry_max=float(other["ry_max"]),
        prior_other=float(other["prior_other"]),
        n_boot=int(boot_cfg.get("n_boot", 2000)),
        seed=int(boot_cfg.get("seed", 0)),
        karyotypes=karyos,
    )


def default_config() -> Config:
    return config_from_mapping(DEFAULT_CONFIG)


def load_config(path: str | Path) -> Config:
    with open(path, "rb") as fh:
        return config_from_mapping(tomllib.load(fh))


def norm(c: str) -> str:
    return c.lower().removeprefix("chr")


def hits(s: int, e: int, r: tuple[int, int]) -> bool:
    return not (e <= r[0] or s >= r[1])


def npdf(x: float, mu: float, sig: float) -> float:
    return math.exp(-0.5 * ((x - mu) / sig) ** 2) / (sig * math.sqrt(2 * math.pi))


def load(path: str, cfg: Config) -> tuple[list[float], list[float], list[float]]:
    """Read mosdepth regions bed(.gz); split depths: autosome / chrX-nonPAR / chrY-euchromatin."""
    auto: list[float] = []
    x: list[float] = []
    y: list[float] = []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c, s, e, d = norm(f[0]), int(f[1]), int(f[2]), float(f[3])
            if c in AUTO:
                auto.append(d)
            elif c == "x":
                if not (hits(s, e, cfg.par1) or hits(s, e, cfg.par2)):
                    x.append(d)
            elif c == "y":
                if hits(s, e, cfg.y_euchromatin):
                    y.append(d)
    return auto, x, y


def boot(xw: list[float], yw: list[float], aw: list[float], cfg: Config):
    """Bootstrap 95% CIs and SEs for Rx and Ry by resampling windows."""
    rxs: list[float] = []
    rys: list[float] = []
    for _ in range(cfg.n_boot):
        a = st.median(random.choices(aw, k=len(aw)))
        if a <= 0:
            continue
        rxs.append(st.median(random.choices(xw, k=len(xw))) / a)
        if yw:
            rys.append(st.median(random.choices(yw, k=len(yw))) / a)
    rxs.sort()
    rys.sort()

    def ci(v: list[float]) -> tuple[float, float]:
        if not v:
            return (float("nan"), float("nan"))
        return (v[int(0.025 * len(v))], v[int(0.975 * len(v))])

    def se(v: list[float]) -> float:
        return st.pstdev(v) if len(v) > 1 else 0.0

    return ci(rxs), se(rxs), ci(rys), se(rys)


def classify(rx: float, ry: float, se_rx: float, se_ry: float, cfg: Config):
    """Posterior over the karyotype set + diffuse OTHER, ranked high→low."""
    sx = math.sqrt(cfg.sigma_rx_model**2 + se_rx**2)
    sy = math.sqrt(cfg.sigma_ry_model**2 + se_ry**2)
    sc: dict[str, float] = {}
    for k in cfg.karyotypes:
        sc[k.label] = (
            k.prior
            * npdf(rx, k.n_x * cfg.rx_unit, sx)
            * npdf(ry, k.n_y * cfg.ry_unit, sy)
        )
    in_box = (0 <= rx <= cfg.rx_max) and (0 <= ry <= cfg.ry_max)
    sc["OTHER"] = cfg.prior_other * (1.0 / (cfg.rx_max * cfg.ry_max) if in_box else 1e-300)
    z = sum(sc.values()) or 1e-300
    post = {k: v / z for k, v in sc.items()}
    return sorted(post.items(), key=lambda t: -t[1])


def rx_ry(auto: list[float], x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Point estimates: autosome median depth, Rx, Ry."""
    a = st.median(auto)
    rx = st.median(x) / a
    ry = (st.median(y) / a) if y else 0.0
    return a, rx, ry


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Karyotype from mosdepth (Rx, Ry) posterior on T2T-CHM13v2.0.",
    )
    ap.add_argument("-c", "--config", help="per-platform TOML config (default: built-in)")
    ap.add_argument("regions", nargs="+", help="mosdepth *.regions.bed(.gz) file(s)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config) if args.config else default_config()
    random.seed(cfg.seed)

    print("sample\tauto_dp\tRx\tRx_CI\tRy\tRy_CI\tcall\tconf\trunner_up")
    for path in args.regions:
        auto, x, y = load(path, cfg)
        if not auto or not x:
            print(f"{path}\tinsufficient windows")
            continue
        a, rx, ry = rx_ry(auto, x, y)
        (rxlo, rxhi), se_rx, (rylo, ryhi), se_ry = boot(x, y, auto, cfg)
        ranked = classify(rx, ry, se_rx, se_ry, cfg)
        (c1, p1), (c2, p2) = ranked[0], ranked[1]
        print(
            f"{path}\t{a:.1f}\t{rx:.3f}\t{rxlo:.3f}-{rxhi:.3f}\t"
            f"{ry:.3f}\t{rylo:.3f}-{ryhi:.3f}\t{c1}\t{p1:.4f}\t{c2}:{p2:.3f}"
        )


if __name__ == "__main__":
    main()
