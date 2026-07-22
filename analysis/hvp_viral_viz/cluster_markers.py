"""Find per-cluster marker ORFs and annotate them with UniProt protein names.

Cluster labels in label_clusters.py come from joins on the taxid (host_group,
ICTV family/order). The clustering signal itself comes from the protein-level
foldseek m8 hits (each ORF query is a 3Di structural-similarity hit to a
BFVD/UniProt protein). This script closes the loop: pull the ORFs that are
characteristic of each cluster, look up what those proteins are.

Inputs:
    <sample>/hits_filtered.parquet
    <sample>/plots/clusters.parquet     written by plots.plot_umap
    refs/cache/uniprot_protein_name.parquet  (optional but expected)

Outputs:
    <sample>/plots/cluster_markers.tsv

For each cluster, picks top-K marker ORFs via scanpy.tl.rank_genes_groups
(Wilcoxon over the virus × ORF count matrix). For each marker ORF, picks the
best (max bits) foldseek hit and joins to protein_name. Default K=10.

Usage:
    PYTHONPATH=src python3 -m hvp_viral_viz.cluster_markers \\
        --sample-dir out/HVP-0006.1_34P
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


_HITS_NAMES = ("hits_filtered.parquet", "pooled_hits.parquet")
_ANNOT_COLS = ("query", "target", "uniprot_id", "bits", "fident", "evalue")


def resolve_hits_path(sample_dir: Path) -> Path:
    """Return the per-sample / pooled hits file in sample_dir, whichever exists.

    Per-sample dirs have hits_filtered.parquet; pool dirs have pooled_hits.parquet.
    Identical schema (the pool is concat of per-sample post-filter frames).
    """
    for name in _HITS_NAMES:
        p = sample_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"no hits parquet in {sample_dir} — tried {' / '.join(_HITS_NAMES)}"
    )


def load_hits_for_matrix(path: Path) -> pd.DataFrame:
    """Load only the (taxid, query) columns needed to build the virus×ORF matrix.

    Column-projection via pyarrow keeps peak memory bounded for pool-scale
    hits files (92M+ rows). Per-sample files are tiny either way.
    """
    import pyarrow.parquet as pq
    return pq.read_table(path, columns=["taxid", "query"]).to_pandas()


def load_hits_filtered_to_queries(path: Path, queries: Iterable[str]) -> pd.DataFrame:
    """Load only rows whose `query` is in `queries`, with annotation columns.

    Used by annotate_markers — typically ~10 markers × N clusters queries
    (a few hundred to a few thousand rows). Pushes the filter into parquet
    so the full DataFrame never materializes.
    """
    import pyarrow.parquet as pq
    q_list = list(queries)
    if not q_list:
        return pd.DataFrame(columns=list(_ANNOT_COLS))
    table = pq.read_table(
        path,
        columns=list(_ANNOT_COLS),
        filters=[("query", "in", q_list)],
    )
    return table.to_pandas()


def compute_marker_orfs(hits: pd.DataFrame, clusters: pd.DataFrame,
                        n_top: int = 10) -> pd.DataFrame:
    """Returns one row per (cluster, marker_rank). Columns: cluster, marker_rank, query, score, pval."""
    import anndata as ad
    import scanpy as sc

    from hvp_viral_viz.plots import _build_virus_orf_matrix

    X, taxids, queries = _build_virus_orf_matrix(hits)
    print(f"[markers] matrix {X.shape[0]:,} viruses × {X.shape[1]:,} ORFs", file=sys.stderr)

    # Build cluster vector aligned to taxids.
    cmap = dict(zip(clusters["taxid"].astype("int64"), clusters["cluster"].astype(str)))
    cl = pd.Series([cmap.get(t, "_unassigned") for t in taxids], dtype=str)
    keep = (cl != "_unassigned").to_numpy()
    if keep.sum() < len(taxids):
        print(f"[markers]   dropping {(~keep).sum()} viruses without cluster assignment",
              file=sys.stderr)

    a = ad.AnnData(
        X=X[keep],
        obs=pd.DataFrame({"leiden": cl[keep].to_numpy()},
                         index=[str(t) for t, k in zip(taxids, keep) if k]),
        var=pd.DataFrame(index=queries),
    )
    a.obs["leiden"] = pd.Categorical(a.obs["leiden"])

    group_sizes = a.obs["leiden"].value_counts()
    drop_groups = group_sizes[group_sizes < 2].index.tolist()
    if drop_groups:
        print(f"[markers]   dropping {len(drop_groups)} singleton clusters "
              f"(rank_genes_groups needs ≥2): {sorted(drop_groups)}",
              file=sys.stderr)
        keep_mask = ~a.obs["leiden"].isin(drop_groups)
        a = a[keep_mask].copy()
        a.obs["leiden"] = pd.Categorical(
            a.obs["leiden"].astype(str),
            categories=sorted(set(a.obs["leiden"].astype(str)))
        )

    print(f"[markers] sc.tl.rank_genes_groups (wilcoxon, top {n_top})", file=sys.stderr)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon", n_genes=n_top,
                            use_raw=False)
    res = a.uns["rank_genes_groups"]

    rows = []
    for cluster in res["names"].dtype.names:
        names = res["names"][cluster]
        scores = res["scores"][cluster]
        pvals = res["pvals_adj"][cluster]
        for rank, (q, s, p) in enumerate(zip(names, scores, pvals)):
            rows.append({
                "cluster": cluster,
                "marker_rank": rank + 1,
                "query": str(q),
                "score": float(s),
                "pval_adj": float(p),
            })
    return pd.DataFrame(rows)


_BAD_NAMES = {"", "deleted", "Deleted", "DELETED", "Uncharacterized protein"}


def annotate_markers(markers: pd.DataFrame, hits: pd.DataFrame,
                     protein_names: pd.DataFrame | None) -> pd.DataFrame:
    """Join each marker query to its best BFVD hit + protein name. When the
    top-bits hit returns "deleted" / no name, walks down hit list to find a
    real annotation."""
    cols = ["query", "target", "uniprot_id", "bits", "fident", "evalue"]
    if protein_names is not None and len(protein_names):
        h = hits[cols].merge(protein_names[["uniprot_id", "protein_name"]],
                             on="uniprot_id", how="left")
    else:
        h = hits[cols].copy()
        h["protein_name"] = ""
    h["protein_name"] = h["protein_name"].fillna("")
    h = h[h["query"].isin(set(markers["query"]))].sort_values(
        ["query", "bits"], ascending=[True, False]
    )
    h["_real"] = ~h["protein_name"].isin(_BAD_NAMES)

    # Top-bits per query.
    top = h.drop_duplicates(subset=["query"], keep="first").copy()
    # Top real-annotated per query (may be empty for some queries).
    real = h[h["_real"]].drop_duplicates(subset=["query"], keep="first")
    real_q = set(real["query"])
    # Merge: use real where available, else fall back to top.
    fallback = top[~top["query"].isin(real_q)]
    best = pd.concat([real, fallback], ignore_index=True).drop(columns=["_real"])

    out = markers.merge(best, on="query", how="left")
    return out


def summarize(annotated: pd.DataFrame, labels: pd.DataFrame | None = None) -> pd.DataFrame:
    """Collapse to one row per cluster with semicolon-joined top markers."""
    summary_rows = []
    for cl, grp in annotated.groupby("cluster"):
        grp = grp.sort_values("marker_rank")
        names = []
        for _, r in grp.iterrows():
            nm = r.get("protein_name") or ""
            if not isinstance(nm, str) or not nm.strip():
                nm = f"(no name) {r.get('uniprot_id', '')}"
            names.append(f"{nm} [b={r.get('bits', 0):.0f}]")
        summary_rows.append({
            "cluster": cl,
            "n_markers": len(grp),
            "top_proteins": " ; ".join(names[:5]),
            "all_uniprots": " ; ".join(str(u) for u in grp["uniprot_id"].dropna()),
        })
    out = pd.DataFrame(summary_rows)
    if labels is not None:
        lab = labels[["cluster", "n_viruses", "label"]].copy()
        lab["cluster"] = lab["cluster"].astype(str)
        out["cluster"] = out["cluster"].astype(str)
        out = lab.merge(out, on="cluster", how="left").sort_values("n_viruses", ascending=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample-dir", type=Path, required=True)
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument("--rank", default="order", choices=["class", "order", "family", "genus"],
                    help="ICTV rank to merge cluster labels from (matches label_clusters.py)")
    ap.add_argument("--n-top", type=int, default=10)
    args = ap.parse_args()

    cl_path = args.sample_dir / "plots" / "clusters.parquet"
    if not cl_path.exists():
        sys.exit(f"[markers] missing {cl_path} — run plots first")
    h_path = resolve_hits_path(args.sample_dir)

    clusters = pd.read_parquet(cl_path)
    hits_matrix = load_hits_for_matrix(h_path)
    print(f"[markers] {len(hits_matrix):,} hits (matrix pass), {len(clusters):,} viruses, "
          f"{clusters['cluster'].nunique()} clusters; source={h_path.name}",
          file=sys.stderr)

    pn_path = args.refs_dir / "cache" / "uniprot_protein_name.parquet"
    if pn_path.exists():
        pn = pd.read_parquet(pn_path)
        print(f"[markers] protein-name cache: {len(pn):,} entries", file=sys.stderr)
    else:
        pn = None
        print(f"[markers] protein-name cache missing → annotation will be UniProt-id only",
              file=sys.stderr)

    markers = compute_marker_orfs(hits_matrix, clusters, n_top=args.n_top)
    del hits_matrix
    hits_annot = load_hits_filtered_to_queries(h_path, set(markers["query"]))
    print(f"[markers] annotation pass: {len(hits_annot):,} hits for "
          f"{markers['query'].nunique()} marker queries", file=sys.stderr)
    annotated = annotate_markers(markers, hits_annot, pn)

    out_tsv = args.sample_dir / "plots" / "cluster_markers.tsv"
    annotated.to_csv(out_tsv, sep="\t", index=False)
    print(f"[markers] wrote {out_tsv}", file=sys.stderr)

    labels_path = args.sample_dir / "plots" / f"cluster_labels_{args.rank}.tsv"
    labels = pd.read_csv(labels_path, sep="\t") if labels_path.exists() else None
    summary = summarize(annotated, labels=labels)
    summary_tsv = args.sample_dir / "plots" / "cluster_markers_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    print(f"[markers] wrote {summary_tsv}", file=sys.stderr)

    # Print markdown summary.
    print()
    print(f"# Per-cluster marker proteins ({len(summary)} clusters)")
    print()
    for _, r in summary.iterrows():
        head = f"## cluster {r['cluster']}"
        if "label" in r and pd.notna(r.get("label")):
            head += f"  —  {r['label']}  (n={int(r.get('n_viruses', 0)):,})"
        print(head)
        print(f"  top: {r['top_proteins']}")
        print()


if __name__ == "__main__":
    main()
