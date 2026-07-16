"""Recolor pooled-cohort plots by specimen type (gut vs oral).

Reads:
- sample_diag_coords.parquet (45 samples × PCA + UMAP coords from cohort_extras)
- virus_samples.parquet (per-virus list of samples it was detected in)
- clusters.parquet (per-virus unified UMAP coords)
- HVP-0006_1.tsv (Terra data table; semibin_environment column =
  human_gut|human_oral — used here as proxy for stool vs saliva)

Writes (into the same plots/ dir):
- sample_diag_pca_by_specimen.png
- sample_diag_umap_by_specimen.png
- umap_per_specimen.png        (3-panel: gut-only / oral-only / shared)
- umap_by_specimen_fraction.png (continuous gut↔oral fraction)
- specimen_labels.tsv           (sample → specimen mapping actually used)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch


SPECIMEN_COLORS = {
    "human_gut": "#d2691e",   # warm brown — stool
    "human_oral": "#1f77b4",  # cool blue — saliva
}
SPECIMEN_LABEL = {
    "human_gut": "stool (human_gut)",
    "human_oral": "saliva (human_oral)",
}


def _save(fig, path: Path, footer: str) -> None:
    fig.text(0.01, 0.005, footer, ha="left", va="bottom",
             fontsize=7, color="gray", family="monospace")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"[specimen]   wrote {path}", file=sys.stderr)


def load_specimen_labels(table_tsv: Path) -> pd.DataFrame:
    """Read HVP-0006_1.tsv → DataFrame(sample, specimen)."""
    df = pd.read_csv(table_tsv, sep="\t", usecols=["name", "semibin_environment"])
    df = df.rename(columns={"name": "sample", "semibin_environment": "specimen"})
    return df


def _scatter_diag(coords: pd.DataFrame, x: str, y: str, title: str,
                  out: Path, footer: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 8.5))
    n_hits = coords["n_hits"].to_numpy()
    size = 40 + 360 * (n_hits - n_hits.min()) / max(1, n_hits.max() - n_hits.min())

    for spec, sub in coords.groupby("specimen"):
        c = SPECIMEN_COLORS.get(spec, "#888888")
        ax.scatter(sub[x], sub[y],
                   s=size[sub.index], c=c, alpha=0.85,
                   edgecolors="black", linewidths=0.4,
                   label=SPECIMEN_LABEL.get(spec, spec))

    for _, r in coords.iterrows():
        label = r["sample"].split("_", 1)[1] if "_" in r["sample"] else r["sample"]
        ax.annotate(label, (r[x], r[y]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", alpha=0.85)

    ax.set_xlabel(x); ax.set_ylabel(y); ax.set_title(title)
    ax.legend(loc="best", frameon=True, fontsize=9)
    ax.text(0.99, 0.01, "marker size ∝ library depth (n_hits)",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, color="#555")
    _save(fig, out, footer)


def plot_sample_diag_by_specimen(coords_pq: Path, specimen: pd.DataFrame,
                                 plots_dir: Path, footer: str) -> None:
    coords = pd.read_parquet(coords_pq).merge(specimen, on="sample", how="left")
    missing = coords["specimen"].isna().sum()
    if missing:
        print(f"[specimen] WARN: {missing} samples have no specimen label",
              file=sys.stderr)
    coords = coords.reset_index(drop=True)

    _scatter_diag(
        coords, "pc1", "pc2",
        f"Cohort structure (PCA) by specimen  —  {len(coords)} samples",
        plots_dir / "sample_diag_pca_by_specimen.png", footer,
    )
    _scatter_diag(
        coords, "umap_x", "umap_y",
        f"Cohort structure (UMAP) by specimen  —  {len(coords)} samples",
        plots_dir / "sample_diag_umap_by_specimen.png", footer,
    )


def annotate_viruses_by_specimen(vs: pd.DataFrame,
                                 sample_to_spec: dict[str, str]) -> pd.DataFrame:
    """For each virus, count #gut samples vs #oral samples."""
    def _counts(samples: np.ndarray) -> tuple[int, int]:
        g = o = 0
        for s in samples:
            spec = sample_to_spec.get(s)
            if spec == "human_gut": g += 1
            elif spec == "human_oral": o += 1
        return g, o

    counts = vs["samples_present"].apply(_counts)
    vs = vs.copy()
    vs["n_gut"] = counts.apply(lambda t: t[0])
    vs["n_oral"] = counts.apply(lambda t: t[1])
    return vs


def plot_umap_per_specimen(clusters: pd.DataFrame, vs: pd.DataFrame,
                           out: Path, footer: str) -> None:
    """3-panel: gut-only, oral-only, shared viruses on the unified UMAP."""
    merged = clusters.merge(vs[["taxid", "n_gut", "n_oral"]],
                            on="taxid", how="left")
    merged[["n_gut", "n_oral"]] = merged[["n_gut", "n_oral"]].fillna(0).astype(int)

    gut_only = (merged["n_gut"] > 0) & (merged["n_oral"] == 0)
    oral_only = (merged["n_oral"] > 0) & (merged["n_gut"] == 0)
    shared = (merged["n_gut"] > 0) & (merged["n_oral"] > 0)
    none = ~(gut_only | oral_only | shared)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    xmin, xmax = merged["umap_x"].min() - 0.5, merged["umap_x"].max() + 0.5
    ymin, ymax = merged["umap_y"].min() - 0.5, merged["umap_y"].max() + 0.5

    panels = [
        ("stool-only", gut_only, SPECIMEN_COLORS["human_gut"]),
        ("saliva-only", oral_only, SPECIMEN_COLORS["human_oral"]),
        ("shared (≥1 of each)", shared, "#6a3d9a"),
    ]
    for ax, (name, mask, color) in zip(axes, panels):
        ax.scatter(merged.loc[~mask, "umap_x"], merged.loc[~mask, "umap_y"],
                   s=2.0, c="#e0e0e0", alpha=0.45, linewidths=0)
        ax.scatter(merged.loc[mask, "umap_x"], merged.loc[mask, "umap_y"],
                   s=3.5, c=color, alpha=0.85, linewidths=0)
        ax.set_title(f"{name}  (n={int(mask.sum()):,})", fontsize=11)
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
        ax.set_xticks([]); ax.set_yticks([])

    n_g, n_o, n_s, n_n = int(gut_only.sum()), int(oral_only.sum()), int(shared.sum()), int(none.sum())
    fig.suptitle(
        f"Unified virus UMAP by specimen-of-origin  —  {len(merged):,} viruses\n"
        f"stool-only={n_g:,}   saliva-only={n_o:,}   shared={n_s:,}   neither={n_n:,}",
        fontsize=12, y=1.005,
    )
    _save(fig, out, footer)


def plot_umap_by_specimen_fraction(clusters: pd.DataFrame, vs: pd.DataFrame,
                                   out: Path, footer: str,
                                   n_gut_samples: int, n_oral_samples: int) -> None:
    """Continuous coloring: rate-normalized gut-fraction per virus.

    Frac = (n_oral/N_oral) / (n_gut/N_gut + n_oral/N_oral)
    → 0 = pure stool, 1 = pure saliva, 0.5 = equal prevalence rate.
    Viruses with zero detections in both groups are dropped.
    """
    merged = clusters.merge(vs[["taxid", "n_gut", "n_oral"]],
                            on="taxid", how="left")
    merged[["n_gut", "n_oral"]] = merged[["n_gut", "n_oral"]].fillna(0).astype(int)

    rate_g = merged["n_gut"] / max(1, n_gut_samples)
    rate_o = merged["n_oral"] / max(1, n_oral_samples)
    denom = rate_g + rate_o
    detected = denom > 0
    frac = pd.Series(np.nan, index=merged.index)
    frac[detected] = rate_o[detected] / denom[detected]

    fig, ax = plt.subplots(figsize=(12, 8.5))
    # Backdrop: undetected (should be empty by construction — all viruses
    # in clusters had ≥1 hit somewhere). Plot them in gray just in case.
    if (~detected).any():
        ax.scatter(merged.loc[~detected, "umap_x"], merged.loc[~detected, "umap_y"],
                   s=2.0, c="#e0e0e0", alpha=0.5, linewidths=0)
    sc = ax.scatter(merged.loc[detected, "umap_x"],
                    merged.loc[detected, "umap_y"],
                    s=4, c=frac[detected], cmap="coolwarm",
                    vmin=0, vmax=1, alpha=0.85, linewidths=0)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Unified virus UMAP — saliva-vs-stool prevalence rate\n"
        f"per-virus rate = n_specimen / N_specimen; "
        f"N_stool={n_gut_samples}, N_saliva={n_oral_samples}"
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.7,
                      label="saliva-fraction  (0 = pure stool · 1 = pure saliva)")
    cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    _save(fig, out, footer)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pooled-dir", required=True, type=Path,
                    help="cohort _pooled/ dir (must contain virus_samples.parquet, "
                         "summary.json, plots/clusters.parquet, "
                         "plots/sample_diag_coords.parquet)")
    ap.add_argument("--table-tsv", required=True, type=Path,
                    help="HVP-0006_1.tsv with semibin_environment column")
    args = ap.parse_args()

    pooled_dir = args.pooled_dir
    plots_dir = pooled_dir / "plots"
    summary = json.loads((pooled_dir / "summary.json").read_text())
    footer = (f"pooled cohort: {summary['n_samples']} samples, "
              f"{summary['n_hits_pooled']:,} hits, "
              f"{summary['n_viruses']:,} viruses")

    spec_df = load_specimen_labels(args.table_tsv)
    print(f"[specimen] table labels: {len(spec_df)} rows, "
          f"specimens={dict(spec_df['specimen'].value_counts())}",
          file=sys.stderr)
    # Restrict to actual cohort samples
    cohort_samples = set(summary["samples"])
    spec_used = spec_df[spec_df["sample"].isin(cohort_samples)].copy()
    spec_used.to_csv(plots_dir / "specimen_labels.tsv", sep="\t", index=False)
    counts = spec_used["specimen"].value_counts().to_dict()
    print(f"[specimen] cohort-restricted labels: {counts}", file=sys.stderr)

    # 1. Sample-level PCA + UMAP recolor
    plot_sample_diag_by_specimen(
        plots_dir / "sample_diag_coords.parquet", spec_used, plots_dir, footer,
    )

    # 2. Virus-level: per-specimen masks + continuous fraction
    clusters = pd.read_parquet(plots_dir / "clusters.parquet")
    vs = pd.read_parquet(pooled_dir / "virus_samples.parquet")
    sample_to_spec = dict(zip(spec_used["sample"], spec_used["specimen"]))
    vs = annotate_viruses_by_specimen(vs, sample_to_spec)

    n_gut = counts.get("human_gut", 0)
    n_oral = counts.get("human_oral", 0)
    plot_umap_per_specimen(clusters, vs, plots_dir / "umap_per_specimen.png", footer)
    plot_umap_by_specimen_fraction(
        clusters, vs, plots_dir / "umap_by_specimen_fraction.png",
        footer, n_gut, n_oral,
    )


if __name__ == "__main__":
    main()
