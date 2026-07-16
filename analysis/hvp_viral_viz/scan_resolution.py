"""Pick a leiden resolution by standardized clustering quality + ICTV alignment.

Problem: resolution=1.0 was arbitrary. Different resolutions land at different
biological hierarchy levels — low resolution merges families into orders /
classes; high resolution splits below genus. We want one resolution whose
partition sits on a SINGLE hierarchy rank, not a mixture.

Approach:
  1. Build embedding + kNN graph once (reuses plots.compute_virus_embedding).
  2. Scan resolutions in {0.1 .. 4.0}, run leiden at each.
  3. For each partition, compute:
       - Intrinsic: modularity, silhouette on SVD embedding.
       - Extrinsic: ARI + NMI vs ICTV class / order / family / genus
         (sklearn.metrics; restricted to viruses where rank is known).
  4. Pick the optimum: the (resolution, rank) pair that maximizes ARI vs
     that single rank, with NMI as tiebreaker. This is the rank-coherent
     partition — clusters correspond to one specific ICTV level.

Outputs:
  out/<sample>/plots/leiden_scan.tsv         per-resolution metric table
  out/<sample>/plots/leiden_scan.png         metric curves vs resolution
  stdout                                      summary + recommendation
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RANKS = ["class", "order", "family", "genus"]
RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0]


def graph_modularity(adata, clusters: np.ndarray) -> float:
    """Newman modularity on the scanpy kNN connectivities graph."""
    try:
        import igraph as ig
        import leidenalg
    except ImportError:
        return float("nan")
    conn = adata.obsp["connectivities"]
    # Build undirected igraph from upper-triangle of connectivities.
    coo = conn.tocoo()
    mask = coo.row < coo.col
    edges = list(zip(coo.row[mask].tolist(), coo.col[mask].tolist()))
    weights = coo.data[mask].tolist()
    g = ig.Graph(n=conn.shape[0], edges=edges, directed=False,
                 edge_attrs={"weight": weights})
    membership = pd.Categorical(clusters).codes.tolist()
    return float(g.modularity(membership, weights="weight"))


def score_partition(adata, clusters: np.ndarray, taxids: list,
                    meta: pd.DataFrame) -> dict:
    """Compute all metrics for one resolution's partition."""
    from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                                 silhouette_score)

    out = {"n_clusters": int(len(set(clusters)))}
    out["modularity"] = graph_modularity(adata, clusters)

    # Silhouette on SVD embedding — sample if too many points (O(n²)).
    Xz = adata.obsm["X_pca"]
    n = Xz.shape[0]
    if out["n_clusters"] >= 2 and out["n_clusters"] < n:
        sample_size = min(2000, n)
        try:
            out["silhouette"] = float(silhouette_score(
                Xz, clusters, metric="cosine",
                sample_size=sample_size, random_state=42,
            ))
        except Exception:
            out["silhouette"] = float("nan")
    else:
        out["silhouette"] = float("nan")

    # ARI / NMI vs each rank — restrict to rows where the rank is annotated.
    by_tx = meta.set_index("taxid")
    cluster_s = pd.Series(clusters, index=taxids)
    for rank in RANKS:
        if rank not in by_tx.columns:
            out[f"ari_{rank}"] = float("nan")
            out[f"nmi_{rank}"] = float("nan")
            out[f"n_{rank}"] = 0
            continue
        rank_s = by_tx[rank].reindex(taxids)
        mask = rank_s.notna() & cluster_s.notna()
        if mask.sum() < 50 or rank_s[mask].nunique() < 2:
            out[f"ari_{rank}"] = float("nan")
            out[f"nmi_{rank}"] = float("nan")
            out[f"n_{rank}"] = int(mask.sum())
            continue
        out[f"ari_{rank}"] = float(adjusted_rand_score(
            rank_s[mask].astype(str), cluster_s[mask].astype(str)
        ))
        out[f"nmi_{rank}"] = float(normalized_mutual_info_score(
            rank_s[mask].astype(str), cluster_s[mask].astype(str)
        ))
        out[f"n_{rank}"] = int(mask.sum())
    return out


def scan(adata, taxids, meta: pd.DataFrame) -> pd.DataFrame:
    import scanpy as sc
    rows = []
    for r in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=r, random_state=42, flavor="igraph",
                     n_iterations=2, directed=False, key_added="_leiden_scan")
        clusters = adata.obs["_leiden_scan"].astype(str).to_numpy()
        scores = score_partition(adata, clusters, taxids, meta)
        scores["resolution"] = r
        rows.append(scores)
        print(f"[scan]   res={r:.2f}  k={scores['n_clusters']:>3}  "
              f"mod={scores['modularity']:.3f}  sil={scores['silhouette']:.3f}  "
              f"ARI(class/order/family/genus)="
              f"{scores['ari_class']:.3f}/{scores['ari_order']:.3f}/"
              f"{scores['ari_family']:.3f}/{scores['ari_genus']:.3f}",
              file=sys.stderr)
    cols = ["resolution", "n_clusters", "modularity", "silhouette"]
    for rank in RANKS:
        cols += [f"ari_{rank}", f"nmi_{rank}", f"n_{rank}"]
    return pd.DataFrame(rows)[cols]


def recommend(table: pd.DataFrame) -> tuple[float, str, dict]:
    """Pick rank-coherent (resolution, rank) by max ARI; NMI tiebreak."""
    best = (float("-inf"), float("-inf"), None, None)
    for _, row in table.iterrows():
        for rank in RANKS:
            ari = row[f"ari_{rank}"]
            nmi = row[f"nmi_{rank}"]
            if not np.isfinite(ari):
                continue
            key = (ari, nmi)
            if key > (best[0], best[1]):
                best = (ari, nmi, row["resolution"], rank)
    return best[2], best[3], {"ari": best[0], "nmi": best[1]}


def plot_curves(table: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax2 = ax.twinx()
    ax.plot(table["resolution"], table["modularity"], "o-", label="modularity", color="#1f77b4")
    ax.plot(table["resolution"], table["silhouette"], "s-", label="silhouette", color="#ff7f0e")
    ax2.plot(table["resolution"], table["n_clusters"], "^--",
             label="n_clusters", color="#888888", alpha=0.7)
    ax.set_xlabel("leiden resolution")
    ax.set_ylabel("intrinsic quality")
    ax2.set_ylabel("n_clusters (right axis)", color="#888888")
    ax.set_xscale("log")
    ax.set_title("Intrinsic cluster quality vs resolution")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    palette = {"class": "#d62728", "order": "#ff7f0e",
               "family": "#2ca02c", "genus": "#1f77b4"}
    for rank in RANKS:
        ax.plot(table["resolution"], table[f"ari_{rank}"], "o-",
                label=f"ARI vs {rank}", color=palette[rank])
    ax.set_xlabel("leiden resolution")
    ax.set_ylabel("Adjusted Rand Index vs ICTV rank")
    ax.set_xscale("log")
    ax.set_title("Extrinsic match to ICTV hierarchy")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Leiden resolution scan — pick rank-coherent partition", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[scan] wrote {out_path}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample-dir", type=Path, required=True)
    args = ap.parse_args()

    import anndata as ad
    from hvp_viral_viz.plots import compute_virus_embedding

    hits = pd.read_parquet(args.sample_dir / "hits_filtered.parquet")
    a = ad.read_h5ad(args.sample_dir / "anndata.h5ad")
    meta = a.var.copy().reset_index(drop=True)
    meta = meta.dropna(subset=["taxid"]).copy()
    meta["taxid"] = meta["taxid"].astype("int64")
    print(f"[scan] {len(hits):,} hits, {len(meta):,} viruses; ranks available: "
          f"{[r for r in RANKS if r in meta.columns]}", file=sys.stderr)

    adata, taxids, _ = compute_virus_embedding(hits)
    print(f"[scan] scanning {len(RESOLUTIONS)} resolutions", file=sys.stderr)
    table = scan(adata, taxids, meta)

    out_tsv = args.sample_dir / "plots" / "leiden_scan.tsv"
    table.to_csv(out_tsv, sep="\t", index=False)
    print(f"[scan] wrote {out_tsv}", file=sys.stderr)
    plot_curves(table, args.sample_dir / "plots" / "leiden_scan.png")

    res, rank, scores = recommend(table)
    print()
    print("=" * 70)
    print(f"RECOMMENDATION: resolution={res}  (rank-coherent at ICTV {rank})")
    print(f"  ARI vs {rank} = {scores['ari']:.3f}")
    print(f"  NMI vs {rank} = {scores['nmi']:.3f}")
    chosen_row = table[table["resolution"] == res].iloc[0]
    print(f"  n_clusters    = {int(chosen_row['n_clusters'])}")
    print(f"  modularity    = {chosen_row['modularity']:.3f}")
    print(f"  silhouette    = {chosen_row['silhouette']:.3f}")
    print("=" * 70)
    print()
    print(table.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
