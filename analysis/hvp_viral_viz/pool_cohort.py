"""Pool per-sample hits into unified cohort viral landscape.

Reads every <out_dir>/<sample>/hits_filtered.parquet under a cohort root,
adds a `sample` column, concatenates to pooled_hits.parquet, then builds
a single unified UMAP where each virus (taxid) is one point and features
are sample-prefixed ORFs (queries).

Why `query`, not BFVD `target`, as the feature: target→taxid is 1:1 by
construction in our pipeline (each BFVD entry carries one taxid; the
lineage join attaches that taxid to every row of that target). So a
target-keyed matrix is essentially the identity — no two viruses share
target columns and the UMAP would have no signal. Queries, however,
fan out across multiple targets/taxids via foldseek's top-k hits — one
ORF hits several BFVD entries from different viruses, creating shared
columns that bridge similar viruses. Across the cohort each sample's
queries form their own columns; the SVD then finds latent factors that
unify recurring cross-virus co-occurrence patterns. Sample-prefix the
query (`<sample>::<query>`) so identical strings from different samples
do not collide.

Memory-aware design: 81M-row pooled frame is too large to hold in pandas
with all 27 columns. After concat we write to disk, then every downstream
step streams the parquet via pyarrow with only the columns it needs.

Outputs under <cohort_dir>/_pooled/:
  - pooled_hits.parquet        union of all hits + `sample` column
  - pooled_anndata.h5ad        obs=samples (N), var=viruses, X=sample×virus counts
  - virus_samples.parquet      taxid → samples_present (list) + n_samples
                               (consumed later for sample-coloring of UMAP)
  - plots/
      - abundance_bar_top30.png
      - abundance_by_family_top30.png
      - heatmap_orf_source_top30.png
      - umap_by_family.png
      - umap_by_host.png
      - umap_by_cluster.png
      - clusters.parquet       virus → umap_x/y/cluster
  - summary.json
"""
from __future__ import annotations
import argparse
import gc
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse

from hvp_viral_viz.plots import (
    plot_abundance_bar,
    plot_abundance_by_family,
)

BATCH_ROWS = 2_000_000

META_COLS = ["scientific_name", "realm", "kingdom", "phylum", "class",
             "order", "family", "genus", "species",
             "host_group", "host_source", "host_confidence"]


def _build_pooled_parquet(sample_dirs: list[Path], pooled_pq: Path) -> None:
    """Concat 45 per-sample parquets to a single pooled parquet on disk.

    Streamed write via pyarrow so peak RAM is bounded by one sample at a time.
    """
    import pyarrow as pa

    writer = None
    total = 0
    t0 = time.time()
    try:
        for d in sample_dirs:
            t = pq.read_table(d / "hits_filtered.parquet")
            sample = d.name
            sample_col = pa.array([sample] * t.num_rows, type=pa.string())
            t = t.append_column("sample", sample_col)
            if writer is None:
                writer = pq.ParquetWriter(pooled_pq, t.schema, compression="snappy")
            writer.write_table(t)
            total += t.num_rows
            print(f"[pool]   {sample}: {t.num_rows:,} hits  (running total {total:,})",
                  file=sys.stderr)
            del t
            gc.collect()
    finally:
        if writer is not None:
            writer.close()
    print(f"[pool] wrote {pooled_pq} ({pooled_pq.stat().st_size/1e9:.2f} GB, "
          f"{total:,} hits, {time.time()-t0:.1f}s)", file=sys.stderr)


def _scan_var_and_membership(pooled_pq: Path, samples: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Streaming pass 1: per-virus count, samples_present, metadata,
    orf_source counts, plus global high-conf count and target cardinality.

    Returns (var, virus_samples, stats).
    """
    cols = ["sample", "taxid", "target", "orf_source", "is_high_conf"] + META_COLS
    per_count: dict[int, int] = defaultdict(int)
    per_samples: dict[int, set[str]] = defaultdict(set)
    per_meta: dict[int, dict] = {}
    per_orf: dict[int, Counter] = defaultdict(Counter)
    target_set: set[str] = set()
    n_hc = 0
    n_total = 0
    pf = pq.ParquetFile(pooled_pq)
    t0 = time.time()
    for i, batch in enumerate(pf.iter_batches(batch_size=BATCH_ROWS, columns=cols)):
        df = batch.to_pandas()
        df = df.dropna(subset=["taxid", "target"])
        df["taxid"] = df["taxid"].astype(np.int64)
        n_total += len(df)
        n_hc += int(df["is_high_conf"].fillna(False).sum())
        target_set.update(df["target"].astype(str).unique().tolist())
        for (tx, src), n in df.groupby(["taxid", "orf_source"]).size().items():
            per_orf[tx][src] += int(n)
        for tx, n in df.groupby("taxid").size().items():
            per_count[int(tx)] += int(n)
        for tx, samps in df.groupby("taxid")["sample"].agg(lambda s: set(s)).items():
            per_samples[int(tx)].update(samps)
        # First-seen metadata per taxid
        firsts = df.drop_duplicates(subset=["taxid"])
        for _, r in firsts.iterrows():
            tx = int(r["taxid"])
            if tx not in per_meta:
                per_meta[tx] = {c: r[c] for c in ["scientific_name"] + META_COLS[1:]}
        if (i + 1) % 5 == 0:
            print(f"[pool]   pass1 batch {i+1}  cum_rows={n_total:,}  "
                  f"({time.time()-t0:.1f}s)", file=sys.stderr)
        del df, batch
    print(f"[pool] pass1: {n_total:,} hits, {len(per_count):,} viruses, "
          f"{len(target_set):,} targets, {n_hc:,} high-conf  "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    rows = []
    for tx, cnt in per_count.items():
        m = per_meta.get(tx, {})
        rows.append({
            "taxid": tx,
            "count": cnt,
            "scientific_name": m.get("scientific_name"),
            **{c: m.get(c) for c in META_COLS[1:]},
        })
    var = pd.DataFrame(rows).sort_values("count", ascending=False).reset_index(drop=True)
    var["taxid_str"] = var["taxid"].astype(str)

    vs_rows = []
    for tx, ss in per_samples.items():
        s_sorted = sorted(ss)
        vs_rows.append({"taxid": tx, "samples_present": s_sorted, "n_samples": len(s_sorted)})
    virus_samples = pd.DataFrame(vs_rows).sort_values("n_samples", ascending=False).reset_index(drop=True)

    # orf_source counts → long DF for heatmap
    orf_rows = []
    for tx, c in per_orf.items():
        for src, n in c.items():
            orf_rows.append({"taxid": tx, "orf_source": src, "n": n})
    orf_long = pd.DataFrame(orf_rows)

    stats = {
        "n_hits": n_total,
        "n_high_conf": n_hc,
        "n_viruses": len(var),
        "n_targets": len(target_set),
        "orf_long": orf_long,
    }
    return var, virus_samples, stats


def _build_sample_virus_matrix(pooled_pq: Path,
                               samples: list[str],
                               taxids: list[int]) -> sparse.csr_matrix:
    """Streaming pass 2: sample × virus hit-count sparse matrix."""
    s_idx = {s: i for i, s in enumerate(samples)}
    t_idx = {int(t): i for i, t in enumerate(taxids)}
    pf = pq.ParquetFile(pooled_pq)
    rows, cols, data = [], [], []
    t0 = time.time()
    for i, batch in enumerate(pf.iter_batches(batch_size=BATCH_ROWS,
                                              columns=["sample", "taxid"])):
        df = batch.to_pandas()
        df = df.dropna(subset=["taxid"])
        df["taxid"] = df["taxid"].astype(np.int64)
        rows.append(df["sample"].map(s_idx).astype(np.int32).to_numpy())
        cols.append(df["taxid"].map(t_idx).astype(np.int32).to_numpy())
        data.append(np.ones(len(df), dtype=np.float32))
        del df, batch
    row = np.concatenate(rows); rows = None
    col = np.concatenate(cols); cols = None
    dat = np.concatenate(data); data = None
    M = sparse.coo_matrix((dat, (row, col)),
                          shape=(len(samples), len(taxids))).tocsr()
    M.sum_duplicates()
    print(f"[pool] sample×virus matrix {M.shape}  nnz={M.nnz:,}  "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)
    return M


def _build_virus_query_matrix(pooled_pq: Path,
                              taxids: list[int]) -> tuple[sparse.csr_matrix, int]:
    """Streaming pass 3: virus × (sample-prefixed) query hit-count matrix.

    Filters mirror per-sample _build_virus_orf_matrix:
      - viruses already discovery-filtered upstream
      - queries must hit ≥ 2 distinct viruses (drop singletons)
    Sample-prefix queries so independent samples produce independent columns.
    """
    t_idx = {int(t): i for i, t in enumerate(taxids)}
    pf = pq.ParquetFile(pooled_pq)

    # Pass A — count distinct viruses per (sample, query) → pick cross-virus columns.
    print("[pool] pass3.A — count distinct viruses per (sample, query)",
          file=sys.stderr)
    q_viruses: dict[str, set[int]] = defaultdict(set)
    t0 = time.time()
    for i, batch in enumerate(pf.iter_batches(batch_size=BATCH_ROWS,
                                              columns=["sample", "taxid", "query"])):
        df = batch.to_pandas().dropna(subset=["taxid", "query"])
        df["taxid"] = df["taxid"].astype(np.int64)
        df["q"] = df["sample"].astype(str) + "::" + df["query"].astype(str)
        for q, tx in df.groupby("q")["taxid"].agg(lambda s: set(s)).items():
            q_viruses[q].update(int(x) for x in tx)
        if (i + 1) % 5 == 0:
            print(f"[pool]   pass3.A batch {i+1}: q_seen={len(q_viruses):,}",
                  file=sys.stderr)
        del df, batch
    keep_qs = [q for q, vs in q_viruses.items() if len(vs) >= 2]
    q_viruses = None
    gc.collect()
    keep_qs.sort()
    g_idx = {q: i for i, q in enumerate(keep_qs)}
    print(f"[pool] kept {len(keep_qs):,} cross-virus queries  "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    # Pass B — build COO using kept columns.
    print("[pool] pass3.B — build virus × query COO", file=sys.stderr)
    rows, cols, data = [], [], []
    t1 = time.time()
    for i, batch in enumerate(pf.iter_batches(batch_size=BATCH_ROWS,
                                              columns=["sample", "taxid", "query"])):
        df = batch.to_pandas().dropna(subset=["taxid", "query"])
        df["taxid"] = df["taxid"].astype(np.int64)
        df["q"] = df["sample"].astype(str) + "::" + df["query"].astype(str)
        df = df[df["q"].isin(g_idx)]
        if df.empty:
            continue
        rows.append(df["taxid"].map(t_idx).astype(np.int32).to_numpy())
        cols.append(df["q"].map(g_idx).astype(np.int32).to_numpy())
        data.append(np.ones(len(df), dtype=np.float32))
        if (i + 1) % 5 == 0:
            print(f"[pool]   pass3.B batch {i+1}", file=sys.stderr)
        del df, batch
    row = np.concatenate(rows); rows = None
    col = np.concatenate(cols); cols = None
    dat = np.concatenate(data); data = None
    X = sparse.coo_matrix((dat, (row, col)),
                          shape=(len(taxids), len(keep_qs))).tocsr()
    X.sum_duplicates()
    print(f"[pool] virus×query matrix {X.shape}  nnz={X.nnz:,}  "
          f"({time.time()-t1:.1f}s)", file=sys.stderr)
    return X, len(keep_qs)


def _embed(X: sparse.csr_matrix, n_components: int = 30, n_neighbors: int = 30,
           random_state: int = 42) -> ad.AnnData:
    """log1p + TruncatedSVD → scanpy neighbors + UMAP."""
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import normalize
    import scanpy as sc

    print(f"[embed] log1p median-normalize on {X.shape}", file=sys.stderr)
    X = X.astype(np.float32)
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    target = np.median(row_sums)
    inv = sparse.diags(target / row_sums)
    Xn = inv @ X
    Xn.data = np.log1p(Xn.data)

    n_comp = min(n_components, min(X.shape) - 1)
    print(f"[embed] TruncatedSVD → {n_comp} components", file=sys.stderr)
    svd = TruncatedSVD(n_components=n_comp, random_state=random_state)
    Xz = svd.fit_transform(Xn)
    Xz = normalize(Xz, norm="l2", axis=1)
    print(f"[embed] explained variance: {100 * svd.explained_variance_ratio_.sum():.1f}%",
          file=sys.stderr)

    adata = ad.AnnData(X=X)
    adata.obsm["X_pca"] = Xz
    nn = min(n_neighbors, max(2, Xz.shape[0] - 1))
    print(f"[embed] scanpy neighbors (n_neighbors={nn}, cosine) + UMAP", file=sys.stderr)
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=nn,
                    metric="cosine", random_state=random_state)
    sc.tl.umap(adata, min_dist=0.1, random_state=random_state)
    return adata


def _plot_heatmap_orf_source(var: pd.DataFrame, orf_long: pd.DataFrame,
                             footer: str, out_path: Path, top_n: int = 30) -> None:
    """Top-N viruses × ORF source heatmap from pre-aggregated orf_long counts."""
    import matplotlib.pyplot as plt

    top = var.nlargest(top_n, "count")
    sources = sorted(orf_long["orf_source"].dropna().unique().tolist())
    pivot = (orf_long[orf_long["taxid"].isin(top["taxid"])]
             .pivot_table(index="taxid", columns="orf_source", values="n", aggfunc="sum")
             .reindex(top["taxid"])
             .fillna(0))
    pivot = pivot.reindex(columns=sources, fill_value=0)
    labels = top["scientific_name"].fillna(top["taxid"].astype(str)).tolist()
    M = np.log1p(pivot.values)
    fig, ax = plt.subplots(figsize=(9, max(6, 0.30 * len(labels))))
    im = ax.imshow(M, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=30, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title(f"Top {top_n} viruses: hits per ORF caller (log scale)")
    fig.colorbar(im, ax=ax, label="hits + 1")
    fig.text(0.01, 0.005, footer, ha="left", va="bottom", fontsize=8,
             color="#555555", wrap=True)
    fig.tight_layout(rect=[0, 0.025, 1, 1])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[plot] wrote {out_path}", file=sys.stderr)


def _plot_umaps(adata: ad.AnnData, taxids: list[int], var: pd.DataFrame,
                X_shape: tuple[int, int], footer: str, plots_dir: Path,
                leiden_resolution: float = 1.0) -> None:
    """Render umap_by_family / umap_by_host / umap_by_cluster + clusters.parquet."""
    import matplotlib.pyplot as plt
    import scanpy as sc

    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=42,
                 flavor="igraph", n_iterations=2, directed=False)
    emb = adata.obsm["X_umap"]
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    n_clust = len(set(clusters))
    print(f"[plot] {n_clust} leiden clusters", file=sys.stderr)

    df = pd.DataFrame({
        "taxid": taxids,
        "umap_x": emb[:, 0],
        "umap_y": emb[:, 1],
        "cluster": clusters,
    })
    rank_cols = [c for c in ["class", "order", "family", "genus"] if c in var.columns]
    meta_cols = ["taxid", "scientific_name", "host_group", "count"] + rank_cols
    meta = var[meta_cols].drop_duplicates(subset=["taxid"])
    df = df.merge(meta, on="taxid", how="left")

    (plots_dir / "clusters.parquet").parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(plots_dir / "clusters.parquet", index=False)
    print(f"[plot] wrote {plots_dir/'clusters.parquet'}", file=sys.stderr)

    from hvp_viral_viz.plots import HOST_PALETTE  # palette only

    def _save(fig, out, footer):
        fig.text(0.01, 0.005, footer, ha="left", va="bottom", fontsize=8,
                 color="#555555", wrap=True)
        fig.tight_layout(rect=[0, 0.025, 1, 1])
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"[plot] wrote {out}", file=sys.stderr)

    # 1. UMAP by family
    df["family"] = df["family"].astype("object")
    top_fams = df.dropna(subset=["family"]).groupby("family")["count"].sum().nlargest(12).index
    df["family_grp"] = df["family"].where(df["family"].isin(top_fams), other="other / unassigned")
    fams = list(top_fams) + ["other / unassigned"]
    palette = plt.cm.tab20.colors
    fig, ax = plt.subplots(figsize=(11, 8))
    for i, fam in enumerate(fams):
        sub = df[df["family_grp"] == fam]
        if len(sub) == 0:
            continue
        ax.scatter(sub["umap_x"], sub["umap_y"], s=14,
                   c=[palette[i % len(palette)]], alpha=0.85,
                   label=f"{fam} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Unified UMAP of {len(df):,} viruses by shared-target profile — colored by family\n"
        f"sparse {X_shape[0]:,}×{X_shape[1]:,} hit matrix, log1p + SVD(30), cosine UMAP"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=9, framealpha=0.85)
    _save(fig, plots_dir / "umap_by_family.png", footer)

    # 2. UMAP by host
    host_order = ["unknown", "archaea", "fungus", "metazoan_other", "protist",
                  "plant", "arthropod", "vertebrate_nonhuman", "bacteria", "human"]
    fig, ax = plt.subplots(figsize=(11, 8))
    df["host_group"] = df["host_group"].fillna("unknown")
    for host in host_order:
        sub = df[df["host_group"] == host]
        if len(sub) == 0:
            continue
        color = HOST_PALETTE.get(host, "#bbbbbb")
        alpha = 0.3 if host == "unknown" else 0.75
        size = 6 if host == "unknown" else 12
        ax.scatter(sub["umap_x"], sub["umap_y"], s=size, c=color, alpha=alpha,
                   label=f"{host} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Unified UMAP of {len(df):,} viruses by shared-target profile — colored by host group\n"
        f"sparse {X_shape[0]:,}×{X_shape[1]:,} hit matrix, log1p + SVD(30), cosine UMAP"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=9, framealpha=0.85)
    _save(fig, plots_dir / "umap_by_host.png", footer)

    # 3. UMAP by cluster
    cluster_sizes = df.groupby("cluster").size().sort_values(ascending=False)
    TOP_K_CLUST = 25
    top_clust = cluster_sizes.head(TOP_K_CLUST).index.tolist()
    fig, ax = plt.subplots(figsize=(11, 8))
    other_mask = ~df["cluster"].isin(top_clust)
    if other_mask.any():
        sub = df[other_mask]
        ax.scatter(sub["umap_x"], sub["umap_y"], s=6, c="#dddddd",
                   alpha=0.35, label=f"other ({len(sub):,})", linewidths=0)
    palette = plt.cm.tab20.colors + plt.cm.tab20b.colors
    for i, c in enumerate(top_clust):
        sub = df[df["cluster"] == c]
        ax.scatter(sub["umap_x"], sub["umap_y"], s=12,
                   c=[palette[i % len(palette)]], alpha=0.85,
                   label=f"cluster {c} (n={len(sub):,})", linewidths=0)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(
        f"Unified UMAP of {len(df):,} viruses — colored by leiden cluster "
        f"(resolution={leiden_resolution:.2f}, k={n_clust})"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=7, framealpha=0.85, ncol=2)
    _save(fig, plots_dir / "umap_by_cluster.png", footer)


def pool(cohort_dir: Path, out_dir: Path, leiden_resolution: float = 1.0) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    sample_dirs = sorted(
        d for d in cohort_dir.iterdir()
        if d.is_dir() and not d.name.startswith("_")
        and (d / "hits_filtered.parquet").exists()
    )
    if not sample_dirs:
        sys.exit(f"[pool] no sample dirs with hits_filtered.parquet under {cohort_dir}")
    samples = [d.name for d in sample_dirs]
    print(f"[pool] {len(samples)} samples", file=sys.stderr)

    pooled_pq = out_dir / "pooled_hits.parquet"
    if not pooled_pq.exists():
        _build_pooled_parquet(sample_dirs, pooled_pq)
    else:
        print(f"[pool] reuse {pooled_pq} ({pooled_pq.stat().st_size/1e9:.2f} GB)",
              file=sys.stderr)

    t_top = time.time()

    var, virus_samples, stats = _scan_var_and_membership(pooled_pq, samples)
    var = var.set_index("taxid_str")
    print(f"[pool] var: {len(var):,} viruses, median samples_present "
          f"= {int(virus_samples.n_samples.median())}", file=sys.stderr)
    vs_pq = out_dir / "virus_samples.parquet"
    virus_samples.to_parquet(vs_pq, index=False)
    print(f"[pool] wrote {vs_pq}", file=sys.stderr)
    orf_long = stats.pop("orf_long")

    # Sample × virus AnnData
    taxids = var["taxid"].astype(int).tolist()
    X_sv = _build_sample_virus_matrix(pooled_pq, samples, taxids)
    obs = pd.DataFrame({"sample": samples}).set_index("sample")
    obs["n_hits"] = np.asarray(X_sv.sum(axis=1)).ravel().astype(np.int64)
    obs["n_viruses"] = np.asarray((X_sv > 0).sum(axis=1)).ravel().astype(np.int64)
    adata_sv = ad.AnnData(X=X_sv, obs=obs, var=var)
    adata_sv.uns["counts"] = {
        "n_raw": int(stats["n_hits"]),
        "n_discovery": int(stats["n_hits"]),
        "n_high_conf": int(stats["n_high_conf"]),
        "n_samples": len(samples),
        "n_viruses": len(var),
        "n_targets": int(stats["n_targets"]),
    }
    adata_sv.uns["samples"] = samples
    out_h5ad = out_dir / "pooled_anndata.h5ad"
    adata_sv.write_h5ad(out_h5ad, compression="gzip")
    print(f"[pool] wrote {out_h5ad} ({out_h5ad.stat().st_size/1e6:.1f} MB)",
          file=sys.stderr)

    # Footer + var-driven plots
    footer = (f"pooled cohort: {len(samples)} samples,  "
              f"{stats['n_hits']:,} hits passing discovery  "
              f"({stats['n_high_conf']:,} high-conf),  "
              f"{len(var):,} viruses,  "
              f"{stats['n_targets']:,} BFVD targets")
    var_for_plot = var.reset_index(drop=False).rename(columns={"taxid_str": "_ix"})
    plot_abundance_bar(var_for_plot, footer, plots_dir / "abundance_bar_top30.png")
    plot_abundance_by_family(var_for_plot, footer, plots_dir / "abundance_by_family_top30.png")
    _plot_heatmap_orf_source(var_for_plot, orf_long, footer,
                             plots_dir / "heatmap_orf_source_top30.png")

    # Unified UMAP — virus × sample-prefixed query
    X_vq, n_q = _build_virus_query_matrix(pooled_pq, taxids)
    rs = np.asarray(X_vq.sum(axis=1)).ravel()
    keep_v = np.where(rs > 0)[0]
    if len(keep_v) < len(taxids):
        print(f"[pool] dropping {len(taxids)-len(keep_v):,} viruses with 0 retained queries",
              file=sys.stderr)
        X_vq = X_vq[keep_v]
        taxids_u = [taxids[i] for i in keep_v]
    else:
        taxids_u = taxids

    adata_emb = _embed(X_vq)
    _plot_umaps(adata_emb, taxids_u, var_for_plot,
                X_shape=X_vq.shape, footer=footer, plots_dir=plots_dir,
                leiden_resolution=leiden_resolution)

    summary = {
        "n_samples": len(samples),
        "n_hits_pooled": int(stats["n_hits"]),
        "n_hits_high_conf": int(stats["n_high_conf"]),
        "n_viruses": int(len(var)),
        "n_viruses_in_umap": int(len(taxids_u)),
        "n_bfvd_targets": int(stats["n_targets"]),
        "n_queries_in_umap": int(X_vq.shape[1]),
        "leiden_resolution": leiden_resolution,
        "samples": samples,
        "host_group_breakdown": (
            var_for_plot.groupby("host_group", dropna=False)["count"].sum().to_dict()
        ),
        "outputs": {
            "pooled_hits": str(pooled_pq.relative_to(out_dir.parent)),
            "pooled_anndata": str(out_h5ad.relative_to(out_dir.parent)),
            "virus_samples": str(vs_pq.relative_to(out_dir.parent)),
        },
    }
    summary["host_group_breakdown"] = {
        (str(k) if k is not None else "unknown"): int(v)
        for k, v in summary["host_group_breakdown"].items()
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"[pool] wrote {out_dir/'summary.json'}", file=sys.stderr)
    print(f"[pool] DONE  total {time.time()-t_top:.1f}s", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cohort-dir", type=Path, required=True,
                    help="Cohort root containing one <sample>/ per sample "
                         "with hits_filtered.parquet")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output dir (default: <cohort-dir>/_pooled)")
    ap.add_argument("--leiden-resolution", type=float, default=1.0)
    args = ap.parse_args()
    out = args.out_dir or args.cohort_dir / "_pooled"
    pool(args.cohort_dir, out, leiden_resolution=args.leiden_resolution)


if __name__ == "__main__":
    main()
