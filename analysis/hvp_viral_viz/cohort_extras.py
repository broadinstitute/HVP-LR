"""Cohort-level extras: sample-structure diagnostic + sample-colored viral UMAPs.

Two analyses driven off the pooled cohort artifacts produced by pool_cohort.py:

1. sample-diag — embed the 45 samples by their virus profiles (sample×virus
   matrix = adata.X in pooled_anndata.h5ad). Yields PCA + UMAP scatters of
   45 points, sized by library depth, colored by virus richness. Answers:
   which samples cluster? Outliers?

2. sample-umap — color the unified virus UMAP (clusters.parquet) by sample
   presence (virus_samples.parquet). Two views:
     a. prevalence — viruses colored by n_samples (continuous viridis).
     b. per-sample — 45-panel small multiples, each highlighting viruses
        present in one sample on top of a gray backdrop.
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

import anndata as ad


def _save(fig, path: Path, footer: str) -> None:
    fig.text(0.01, 0.005, footer, ha="left", va="bottom",
             fontsize=7, color="gray", family="monospace")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"[extras]   wrote {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# 1. sample-diag — 45-point cohort structure embedding
# ---------------------------------------------------------------------------

def _embed_samples(X, n_components: int = 10, n_neighbors: int = 10,
                   random_state: int = 42):
    """Log1p + median-normalize sample×virus → PCA + UMAP.

    Returns (Xz_pca, Xz_umap, explained_var).
    """
    from scipy import sparse
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import normalize

    X = X.astype(np.float32)
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    target = np.median(row_sums)
    inv = sparse.diags(target / row_sums)
    Xn = inv @ X
    Xn.data = np.log1p(Xn.data)

    n_comp = min(n_components, min(X.shape) - 1)
    print(f"[extras] TruncatedSVD → {n_comp} components", file=sys.stderr)
    svd = TruncatedSVD(n_components=n_comp, random_state=random_state)
    Xz = svd.fit_transform(Xn)
    Xz = normalize(Xz, norm="l2", axis=1)
    explained = float(svd.explained_variance_ratio_.sum())
    print(f"[extras]   explained variance: {100 * explained:.1f}%", file=sys.stderr)

    import umap
    nn = min(n_neighbors, X.shape[0] - 1)
    print(f"[extras] UMAP (n_neighbors={nn}, cosine, min_dist=0.3)", file=sys.stderr)
    reducer = umap.UMAP(n_neighbors=nn, min_dist=0.3, metric="cosine",
                        random_state=random_state)
    Xu = reducer.fit_transform(Xz)
    return Xz, Xu, explained


def _hierarchical_order(Xz, sample_names):
    """Hierarchical clustering on PCA-reduced sample reps. Returns ordered list of names."""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    d = pdist(Xz, metric="cosine")
    Z = linkage(d, method="average")
    order = leaves_list(Z)
    return [sample_names[i] for i in order], Z


def _plot_sample_diag(adata: ad.AnnData, plots_dir: Path, footer: str) -> dict:
    """PCA scatter + UMAP scatter + dendrogram of 45 samples."""
    samples = list(adata.obs_names)
    n_hits = adata.obs["n_hits"].to_numpy()
    n_vir = adata.obs["n_viruses"].to_numpy()

    Xz, Xu, explained = _embed_samples(adata.X)

    # Marker size scaled by library depth (n_hits), color by n_viruses.
    size_scale = 40 + 360 * (n_hits - n_hits.min()) / max(1, n_hits.max() - n_hits.min())
    norm = Normalize(vmin=n_vir.min(), vmax=n_vir.max())
    cmap = plt.cm.viridis

    # --- PC1-PC2 scatter --------------------------------------------------
    fig, ax = plt.subplots(figsize=(11, 8.5))
    sc = ax.scatter(Xz[:, 0], Xz[:, 1], s=size_scale, c=n_vir, cmap=cmap,
                    norm=norm, alpha=0.85, edgecolors="black", linewidths=0.4)
    for i, s in enumerate(samples):
        # Strip cohort prefix for readability.
        label = s.split("_", 1)[1] if "_" in s else s
        ax.annotate(label, (Xz[i, 0], Xz[i, 1]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", alpha=0.85)
    ax.set_xlabel("PC1 (cosine-normalized)")
    ax.set_ylabel("PC2")
    ax.set_title(
        f"Cohort structure: {len(samples)} samples by virus profile  —  PCA\n"
        f"sample × virus log1p hits, SVD(10) — explained var {100*explained:.1f}%"
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.7, label="distinct viruses per sample")
    ax.text(0.99, 0.01, "marker size ∝ library depth (n_hits)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#555")
    _save(fig, plots_dir / "sample_diag_pca.png", footer)

    # --- UMAP scatter ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(11, 8.5))
    sc = ax.scatter(Xu[:, 0], Xu[:, 1], s=size_scale, c=n_vir, cmap=cmap,
                    norm=norm, alpha=0.85, edgecolors="black", linewidths=0.4)
    for i, s in enumerate(samples):
        label = s.split("_", 1)[1] if "_" in s else s
        ax.annotate(label, (Xu[i, 0], Xu[i, 1]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", alpha=0.85)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Cohort structure: {len(samples)} samples by virus profile  —  UMAP\n"
        f"SVD(10) → cosine UMAP (n_neighbors=10, min_dist=0.3)"
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.7, label="distinct viruses per sample")
    ax.text(0.99, 0.01, "marker size ∝ library depth (n_hits)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#555")
    _save(fig, plots_dir / "sample_diag_umap.png", footer)

    # --- Dendrogram --------------------------------------------------------
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist
    d = pdist(Xz, metric="cosine")
    Z = linkage(d, method="average")
    fig, ax = plt.subplots(figsize=(14, 6))
    short = [s.split("_", 1)[1] if "_" in s else s for s in samples]
    dendrogram(Z, labels=short, leaf_rotation=90, leaf_font_size=8, ax=ax,
               color_threshold=0.6 * max(Z[:, 2]))
    ax.set_ylabel("cosine distance (avg linkage on SVD(10) coords)")
    ax.set_title(f"Cohort dendrogram — {len(samples)} samples by virus profile")
    _save(fig, plots_dir / "sample_diag_dendrogram.png", footer)

    # --- Persist coords ----------------------------------------------------
    df = pd.DataFrame({
        "sample": samples,
        "n_hits": n_hits,
        "n_viruses": n_vir,
        "pc1": Xz[:, 0],
        "pc2": Xz[:, 1],
        "umap_x": Xu[:, 0],
        "umap_y": Xu[:, 1],
    })
    df.to_parquet(plots_dir / "sample_diag_coords.parquet", index=False)
    print(f"[extras]   wrote {plots_dir / 'sample_diag_coords.parquet'}", file=sys.stderr)

    return {"explained_variance": explained}


def main_sample_diag(pooled_dir: Path) -> None:
    h5ad = pooled_dir / "pooled_anndata.h5ad"
    summary_path = pooled_dir / "summary.json"
    plots_dir = pooled_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"[extras] loading {h5ad}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad)
    print(f"[extras]   X: {adata.X.shape} {type(adata.X).__name__}", file=sys.stderr)

    summary = json.loads(summary_path.read_text())
    footer = (f"pooled cohort: {summary['n_samples']} samples, "
              f"{summary['n_hits_pooled']:,} hits, "
              f"{summary['n_viruses']:,} viruses")
    _plot_sample_diag(adata, plots_dir, footer)


# ---------------------------------------------------------------------------
# 2. sample-umap — sample-colored unified viral UMAP
# ---------------------------------------------------------------------------

def _plot_prevalence(clusters: pd.DataFrame, footer: str, out: Path) -> None:
    """Single panel: virus UMAP colored by n_samples present (1..45)."""
    fig, ax = plt.subplots(figsize=(11, 8.5))
    n = clusters["n_samples"].to_numpy()
    sc = ax.scatter(clusters["umap_x"], clusters["umap_y"], s=8,
                    c=n, cmap="viridis", alpha=0.75, linewidths=0)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Unified UMAP of {len(clusters):,} viruses  —  colored by prevalence\n"
        f"per-virus sample presence across cohort (min=1, max={int(n.max())})"
    )
    fig.colorbar(sc, ax=ax, shrink=0.7, label="samples in which virus is detected")
    _save(fig, out, footer)


def _plot_per_sample(clusters: pd.DataFrame, virus_samples: pd.DataFrame,
                     samples: list[str], footer: str, out: Path) -> None:
    """45-panel small multiples: each panel highlights one sample's viruses.

    Backdrop: all 13,822 mapped viruses in light gray.
    Highlight: viruses with the sample in samples_present (colored).
    """
    # Build sample → set of taxids for O(1) membership.
    sample_to_tx: dict[str, set[int]] = {s: set() for s in samples}
    for tx, ss in zip(virus_samples["taxid"], virus_samples["samples_present"]):
        for s in ss:
            if s in sample_to_tx:
                sample_to_tx[s].add(int(tx))

    n_samples = len(samples)
    n_cols = 9
    n_rows = (n_samples + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2.6, n_rows * 2.4),
                             sharex=True, sharey=True)
    axes = axes.flatten()

    palette = (plt.cm.tab20.colors + plt.cm.tab20b.colors
               + plt.cm.tab20c.colors)
    xmin = clusters["umap_x"].min() - 0.5
    xmax = clusters["umap_x"].max() + 0.5
    ymin = clusters["umap_y"].min() - 0.5
    ymax = clusters["umap_y"].max() + 0.5

    for i, s in enumerate(samples):
        ax = axes[i]
        present = clusters["taxid"].isin(sample_to_tx[s])
        absent = ~present
        # Backdrop
        ax.scatter(clusters.loc[absent, "umap_x"],
                   clusters.loc[absent, "umap_y"],
                   s=1.2, c="#e0e0e0", alpha=0.45, linewidths=0)
        # Highlight
        ax.scatter(clusters.loc[present, "umap_x"],
                   clusters.loc[present, "umap_y"],
                   s=2.5, c=[palette[i % len(palette)]], alpha=0.85,
                   linewidths=0)
        label = s.split("_", 1)[1] if "_" in s else s
        ax.set_title(f"{label}  (n={int(present.sum()):,})", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
        for spine in ax.spines.values():
            spine.set_color("#bbbbbb")
            spine.set_linewidth(0.5)

    # Hide unused axes.
    for j in range(n_samples, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"Unified virus UMAP per sample  —  {n_samples} panels  ({len(clusters):,} viruses each)\n"
        "colored = virus detected in that sample, gray = absent",
        fontsize=11, y=0.995,
    )
    _save(fig, out, footer)


def main_sample_umap(pooled_dir: Path) -> None:
    plots_dir = pooled_dir / "plots"
    clusters_pq = plots_dir / "clusters.parquet"
    vs_pq = pooled_dir / "virus_samples.parquet"
    summary_path = pooled_dir / "summary.json"

    print(f"[extras] loading {clusters_pq}", file=sys.stderr)
    clusters = pd.read_parquet(clusters_pq)
    print(f"[extras]   {len(clusters):,} viruses with UMAP coords", file=sys.stderr)

    print(f"[extras] loading {vs_pq}", file=sys.stderr)
    vs = pd.read_parquet(vs_pq)
    # Join n_samples onto clusters.
    clusters = clusters.merge(vs[["taxid", "n_samples"]], on="taxid", how="left")
    clusters["n_samples"] = clusters["n_samples"].fillna(0).astype(int)

    summary = json.loads(summary_path.read_text())
    samples = sorted(summary["samples"])
    footer = (f"pooled cohort: {summary['n_samples']} samples, "
              f"{summary['n_hits_pooled']:,} hits, "
              f"{summary['n_viruses']:,} viruses")

    _plot_prevalence(clusters, footer, plots_dir / "umap_by_prevalence.png")
    _plot_per_sample(clusters, vs, samples, footer,
                     plots_dir / "umap_per_sample.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pooled-dir", required=True, type=Path,
                    help="Path to cohort _pooled/ directory containing pooled_anndata.h5ad, "
                         "virus_samples.parquet, summary.json, and plots/clusters.parquet")
    ap.add_argument("mode", choices=["sample-diag", "sample-umap", "all"],
                    help="sample-diag: 45-point cohort structure UMAP/PCA. "
                         "sample-umap: sample-colored unified virus UMAP "
                         "(prevalence + per-sample small multiples). "
                         "all: run both.")
    args = ap.parse_args()

    if args.mode in ("sample-diag", "all"):
        main_sample_diag(args.pooled_dir)
    if args.mode in ("sample-umap", "all"):
        main_sample_umap(args.pooled_dir)


if __name__ == "__main__":
    main()
