"""Tier 1 visualizations for a single ingested sample.

Reads <out_dir>/<sample>/anndata.h5ad and hits_filtered.parquet,
produces a set of PNG + interactive HTML figures under <out_dir>/<sample>/plots/.

Plots:
    1. abundance_bar           top-N viruses by hit count, colored by host_group
    2. abundance_by_family     hit count summed per family, colored by host_group
    3. heatmap_orf_source      top-N viruses × ORF source (vs2 / assembly / genomad / rescued)
    4. umap_taxonomy           UMAP over per-virus hit-stats, colored by family
    5. umap_host               same UMAP, colored by host_group
    6. footer for every plot:  filter line + N hits passing
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import anndata as ad

TOP_N_BAR = 30
TOP_N_HEATMAP = 30
MIN_HITS_FOR_UMAP = 3  # drop singleton viruses from UMAP (noise)

# Host group palette (consistent across all plots)
HOST_PALETTE = {
    "human": "#d62728",
    "vertebrate_nonhuman": "#ff7f0e",
    "arthropod": "#9467bd",
    "metazoan_other": "#8c564b",
    "plant": "#2ca02c",
    "fungus": "#bcbd22",
    "protist": "#17becf",
    "bacteria": "#1f77b4",
    "archaea": "#7f7f7f",
    "unknown": "#cccccc",
}


def filter_footer(adata: ad.AnnData) -> str:
    thr = adata.uns["thresholds"]["discovery"]
    counts = adata.uns["counts"]
    return (
        f"filter: evalue ≤ {thr['evalue_max']:.0e} AND bits ≥ {thr['bits_min']} "
        f"AND alnlen ≥ {thr['alnlen_min']}    "
        f"N hits passing: {counts['n_discovery']:,} / {counts['n_raw']:,}"
    )


def _save(fig, path: Path, footer: str) -> None:
    fig.text(0.01, 0.005, footer, ha="left", va="bottom",
             fontsize=7, color="gray", family="monospace")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"[plot]   wrote {path}", file=sys.stderr)


def plot_abundance_bar(var: pd.DataFrame, footer: str, out: Path, top_n: int = TOP_N_BAR) -> None:
    df = var.nlargest(top_n, "count").iloc[::-1]
    labels = df.apply(
        lambda r: f"{r['scientific_name']}" if pd.notna(r['scientific_name']) else f"taxid {r['taxid']}",
        axis=1,
    ).str.slice(0, 60)
    colors = df["host_group"].map(HOST_PALETTE).fillna("#cccccc")

    fig, ax = plt.subplots(figsize=(10, max(6, 0.25 * top_n)))
    ax.barh(range(len(df)), df["count"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("foldseek hits passing discovery threshold")
    ax.set_title(f"Top {top_n} viruses by structural hits — single sample")

    # Legend
    used = sorted(set(df["host_group"].dropna()))
    handles = [plt.Rectangle((0, 0), 1, 1, color=HOST_PALETTE.get(h, "#cccccc")) for h in used]
    ax.legend(handles, used, loc="lower right", fontsize=8, title="host group")
    _save(fig, out, footer)


def plot_abundance_by_family(var: pd.DataFrame, footer: str, out: Path, top_n: int = TOP_N_BAR) -> None:
    # Sum counts per family, keep top N families
    fam = (
        var.dropna(subset=["family"])
        .groupby("family")
        .agg(count=("count", "sum"),
             host_group=("host_group", lambda s: s.value_counts().index[0]))
        .nlargest(top_n, "count")
        .iloc[::-1]
    )
    colors = fam["host_group"].map(HOST_PALETTE).fillna("#cccccc")
    fig, ax = plt.subplots(figsize=(10, max(6, 0.25 * top_n)))
    ax.barh(range(len(fam)), fam["count"], color=colors)
    ax.set_yticks(range(len(fam)))
    ax.set_yticklabels(fam.index, fontsize=9)
    ax.set_xlabel("foldseek hits passing discovery threshold")
    ax.set_title(f"Top {top_n} virus families by structural hits — single sample")

    used = sorted(set(fam["host_group"].dropna()))
    handles = [plt.Rectangle((0, 0), 1, 1, color=HOST_PALETTE.get(h, "#cccccc")) for h in used]
    ax.legend(handles, used, loc="lower right", fontsize=8, title="majority host")
    _save(fig, out, footer)


def plot_heatmap_orf_source(hits: pd.DataFrame, var: pd.DataFrame, footer: str,
                             out: Path, top_n: int = TOP_N_HEATMAP) -> None:
    top = var.nlargest(top_n, "count")
    sub = hits[hits["taxid"].isin(top["taxid"])].copy()
    sub["taxid_str"] = sub["taxid"].astype(str)
    mat = sub.groupby(["taxid_str", "orf_source"]).size().unstack(fill_value=0)

    # Order rows by total hits desc
    ordering = top.set_index(top["taxid"].astype(str)).index
    mat = mat.reindex(ordering, fill_value=0)
    # Show all orf sources present
    src_order = [c for c in ["assembly", "vs2", "genomad", "rescued", "unknown"] if c in mat.columns]
    mat = mat[src_order]

    labels = [
        f"{top.loc[top['taxid'].astype(str) == t, 'scientific_name'].iloc[0][:50]}"
        if pd.notna(top.loc[top['taxid'].astype(str) == t, 'scientific_name'].iloc[0])
        else f"taxid {t}"
        for t in ordering
    ]

    fig, ax = plt.subplots(figsize=(6, max(6, 0.3 * top_n)))
    im = ax.imshow(mat.values + 1, aspect="auto", cmap="viridis", norm=LogNorm())
    ax.set_xticks(range(len(src_order)))
    ax.set_xticklabels(src_order, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(f"Top {top_n} viruses: hits per ORF caller (log scale)")
    fig.colorbar(im, ax=ax, label="hits + 1")
    _save(fig, out, footer)


def _build_virus_orf_matrix(hits: pd.DataFrame) -> tuple["sparse.csr_matrix", list[int], list[str]]:
    """Sparse virus × ORF count matrix.

    Each row is a viral taxid, each column is a query ORF. Cells count hits.
    Filters:
      - viruses with ≥ MIN_HITS_FOR_UMAP hits (drops singleton noise)
      - ORFs hitting ≥ 2 distinct viruses (cross-virus signal — singletons
        only add noise to UMAP since they contribute one row each)

    This is the "single-cell-style" representation: viruses ≈ cells,
    ORFs ≈ genes, hit counts ≈ expression. Two viruses are similar when
    their ORF coverage profiles overlap.
    """
    from scipy import sparse

    v_counts = hits.groupby("taxid").size()
    keep_v = v_counts[v_counts >= MIN_HITS_FOR_UMAP].index
    sub = hits[hits["taxid"].isin(keep_v)]

    q_to_v = sub.groupby("query")["taxid"].nunique()
    keep_q = q_to_v[q_to_v >= 2].index
    sub = sub[sub["query"].isin(keep_q)]

    taxids = sorted(sub["taxid"].unique())
    queries = sorted(sub["query"].unique())
    t_idx = {t: i for i, t in enumerate(taxids)}
    q_idx = {q: i for i, q in enumerate(queries)}

    row = sub["taxid"].map(t_idx).values
    col = sub["query"].map(q_idx).values
    data = np.ones(len(sub), dtype=np.float32)
    X = sparse.coo_matrix(
        (data, (row, col)), shape=(len(taxids), len(queries))
    ).tocsr()
    # Collapse duplicate (taxid, query) pairs to count
    X.sum_duplicates()
    return X, taxids, queries


def compute_virus_embedding(hits: pd.DataFrame, n_neighbors: int = 30,
                            n_components: int = 30, random_state: int = 42,
                            verbose: bool = True):
    """Build virus×ORF matrix → normalize → SVD → kNN graph → UMAP.

    Returns (adata, taxids, queries). adata has obsm["X_pca"], obsm["X_umap"],
    and a populated neighbors graph. Does NOT call leiden — caller picks
    resolution.
    """
    import anndata as _ad
    import scanpy as sc
    from scipy import sparse
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import normalize

    if verbose:
        print(f"[embed] building virus × ORF sparse matrix", file=sys.stderr)
    X, taxids, queries = _build_virus_orf_matrix(hits)
    nnz = X.nnz
    if verbose:
        print(
            f"[embed]   shape={X.shape}  nnz={nnz:,}  density={100 * nnz / (X.shape[0] * X.shape[1]):.3f}%",
            file=sys.stderr,
        )
    X = X.astype(np.float32)
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    target = np.median(row_sums)
    inv = sparse.diags(target / row_sums)
    Xn = inv @ X
    Xn.data = np.log1p(Xn.data)

    n_comp = min(n_components, min(X.shape) - 1)
    if verbose:
        print(f"[embed] TruncatedSVD → {n_comp} components", file=sys.stderr)
    svd = TruncatedSVD(n_components=n_comp, random_state=random_state)
    Xz = svd.fit_transform(Xn)
    Xz = normalize(Xz, norm="l2", axis=1)
    if verbose:
        print(f"[embed]   explained variance: {100 * svd.explained_variance_ratio_.sum():.1f}%",
              file=sys.stderr)

    nn = min(n_neighbors, max(2, Xz.shape[0] - 1))
    if verbose:
        print(f"[embed] scanpy neighbors (n_neighbors={nn}, metric=cosine) + UMAP",
              file=sys.stderr)
    adata = _ad.AnnData(X=Xn, obs=pd.DataFrame(index=[str(t) for t in taxids]))
    adata.obsm["X_pca"] = Xz
    sc.pp.neighbors(adata, n_neighbors=nn, use_rep="X_pca",
                    metric="cosine", random_state=random_state)
    sc.tl.umap(adata, min_dist=0.1, random_state=random_state)
    return adata, taxids, queries


def plot_umap(hits: pd.DataFrame, var: pd.DataFrame, footer: str,
              umap_taxonomy_out: Path, umap_host_out: Path,
              umap_cluster_out: Path, leiden_resolution: float = 1.0) -> None:
    import scanpy as sc

    adata, taxids, queries = compute_virus_embedding(hits)
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=42,
                 flavor="igraph", n_iterations=2, directed=False)
    emb = adata.obsm["X_umap"]
    X = adata.X
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    n_clust = len(set(clusters))
    print(f"[plot]   leiden(res={leiden_resolution}): {n_clust} clusters",
          file=sys.stderr)

    key = pd.DataFrame({
        "taxid": taxids,
        "umap_x": emb[:, 0],
        "umap_y": emb[:, 1],
        "cluster": clusters,
    })
    rank_cols = [c for c in ["class", "order", "family", "genus"] if c in var.columns]
    meta_cols = ["taxid", "scientific_name", "host_group", "count"] + rank_cols
    meta = var[meta_cols].drop_duplicates(subset=["taxid"])
    df = key.merge(meta, on="taxid", how="left")

    # Persist cluster assignments for downstream labeling / inspection.
    clusters_out = umap_taxonomy_out.parent / "clusters.parquet"
    df.to_parquet(clusters_out, index=False)
    print(f"[plot]   wrote {clusters_out}", file=sys.stderr)

    # Plot 1 — colored by family (top 12 by hit count, rest grouped 'other')
    df["family"] = df["family"].astype("object")
    top_fams = df.dropna(subset=["family"]).groupby("family")["count"].sum().nlargest(12).index
    df["family_grp"] = df["family"].where(df["family"].isin(top_fams), other="other")

    fig, ax = plt.subplots(figsize=(11, 8.5))
    other = df[df["family_grp"] == "other"]
    ax.scatter(other["umap_x"], other["umap_y"], s=6, c="#dddddd",
               alpha=0.35, label=f"other (n={len(other):,})", linewidths=0)
    palette = plt.cm.tab20.colors
    for i, fam in enumerate(top_fams):
        sub = df[df["family_grp"] == fam]
        ax.scatter(sub["umap_x"], sub["umap_y"], s=14,
                   c=[palette[i % len(palette)]], alpha=0.85,
                   label=f"{fam} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"UMAP of {len(df):,} viruses by shared-ORF profile  —  colored by family\n"
        f"sparse {X.shape[0]:,}×{X.shape[1]:,} hit matrix, log1p + SVD(30), cosine UMAP"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=8, ncol=1, framealpha=0.85)
    _save(fig, umap_taxonomy_out, footer)

    # Plot 2 — colored by host_group; plot 'unknown' first so colored hosts overlay
    fig, ax = plt.subplots(figsize=(11, 8.5))
    host_order = ["unknown", "archaea", "fungus", "metazoan_other", "protist",
                  "plant", "arthropod", "vertebrate_nonhuman", "bacteria", "human"]
    for host in host_order:
        color = HOST_PALETTE.get(host, "#cccccc")
        sub = df[df["host_group"] == host]
        if len(sub) == 0:
            continue
        alpha = 0.3 if host == "unknown" else 0.75
        size = 6 if host == "unknown" else 12
        ax.scatter(sub["umap_x"], sub["umap_y"], s=size, c=color, alpha=alpha,
                   label=f"{host} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"UMAP of {len(df):,} viruses by shared-ORF profile  —  colored by host group\n"
        f"sparse {X.shape[0]:,}×{X.shape[1]:,} hit matrix, log1p + SVD(30), cosine UMAP"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=9, framealpha=0.85)
    _save(fig, umap_host_out, footer)

    # Plot 3 — colored by leiden cluster. Cluster id is unsupervised; sized by
    # cluster member count so big clusters dominate the legend visually.
    fig, ax = plt.subplots(figsize=(11, 8.5))
    cluster_sizes = df.groupby("cluster").size().sort_values(ascending=False)
    palette = plt.cm.tab20.colors + plt.cm.tab20b.colors + plt.cm.tab20c.colors
    # Top-K clusters get distinct colors; small tail merged to "other".
    TOP_K_CLUST = 20
    top_clust = cluster_sizes.head(TOP_K_CLUST).index.tolist()
    other_mask = ~df["cluster"].isin(top_clust)
    if other_mask.any():
        sub = df[other_mask]
        ax.scatter(sub["umap_x"], sub["umap_y"], s=6, c="#dddddd",
                   alpha=0.35, label=f"other ({len(top_clust)}+ tail, n={len(sub):,})",
                   linewidths=0)
    for i, c in enumerate(top_clust):
        sub = df[df["cluster"] == c]
        ax.scatter(sub["umap_x"], sub["umap_y"], s=14,
                   c=[palette[i % len(palette)]], alpha=0.85,
                   label=f"c{c} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"UMAP of {len(df):,} viruses by shared-ORF profile  —  leiden clusters\n"
        f"{n_clust} clusters @ resolution={leiden_resolution} on shared kNN graph (cosine)"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=7, ncol=1, framealpha=0.85)
    _save(fig, umap_cluster_out, footer)


def make_all_plots(sample_dir: Path) -> None:
    h5ad = sample_dir / "anndata.h5ad"
    hits_pq = sample_dir / "hits_filtered.parquet"
    plots_dir = sample_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print(f"[plot] loading {h5ad}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad)
    var = adata.var.reset_index().rename(columns={"index": "taxid_str"})
    print(f"[plot]   {len(var):,} viruses", file=sys.stderr)

    print(f"[plot] loading {hits_pq}", file=sys.stderr)
    hits = pd.read_parquet(hits_pq)
    print(f"[plot]   {len(hits):,} hits", file=sys.stderr)

    footer = filter_footer(adata)

    plot_abundance_bar(var, footer, plots_dir / "abundance_bar_top30.png")
    plot_abundance_by_family(var, footer, plots_dir / "abundance_by_family_top30.png")
    plot_heatmap_orf_source(hits, var, footer, plots_dir / "heatmap_orf_source_top30.png")
    plot_umap(
        hits, var, footer,
        umap_taxonomy_out=plots_dir / "umap_by_family.png",
        umap_host_out=plots_dir / "umap_by_host.png",
        umap_cluster_out=plots_dir / "umap_by_cluster.png",
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample-dir", type=Path, required=True,
                    help="Directory containing anndata.h5ad + hits_filtered.parquet")
    args = ap.parse_args()
    make_all_plots(args.sample_dir)


if __name__ == "__main__":
    main()
