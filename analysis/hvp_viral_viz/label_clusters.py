"""Assign biologically interpretable labels to leiden clusters.

Reads `plots/clusters.parquet` (written by plot_umap) which carries one row
per virus with taxid, cluster id, umap coords, family, host_group,
scientific_name, and the hit count from the sample.

For each cluster, computes:
  - cluster size (n viruses)
  - dominant family + purity (fraction of viruses in that family)
  - dominant host_group + purity
  - 3 representative scientific names (highest hit-count, breaks ties
    by alpha order)
  - secondary family/host if dominant <70%

Emits:
  - cluster_labels.tsv    one row per cluster, sorted by size
  - prints a markdown-ish summary to stdout

Usage:
  PYTHONPATH=src python3 -m hvp_viral_viz.label_clusters \
      --sample-dir out/HVP-0006.1_34P
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd


def _top_label(s: pd.Series, n: int = 2) -> list[tuple[str, float]]:
    """Top-n labels with their fractional share, drop NaN."""
    vc = s.dropna().value_counts(normalize=True)
    return [(str(k), float(v)) for k, v in vc.head(n).items()]


def _format_share(pairs: list[tuple[str, float]]) -> str:
    return ", ".join(f"{k} ({100 * v:.0f}%)" for k, v in pairs)


def _suggest_label_with_rank(rank_top: list[tuple[str, float]],
                             host_top: list[tuple[str, float]],
                             rank_name: str) -> str:
    return _suggest_label(rank_top, host_top)


def _suggest_label(family_top: list[tuple[str, float]],
                   host_top: list[tuple[str, float]]) -> str:
    fam, fam_p = (family_top[0] if family_top else ("unknown", 0.0))
    host, host_p = (host_top[0] if host_top else ("unknown", 0.0))
    # Family qualifiers reflect actual purity, not a fixed 50% cliff.
    if fam_p >= 0.5:
        fam_label = fam
    elif fam_p >= 0.25:
        fam_label = f"{fam}-leaning"
    else:
        fam_label = "polyphyletic"
    if host_p >= 0.7:
        host_label = host
    elif host_p >= 0.4:
        host_label = f"{host}-leaning"
    else:
        host_label = "mixed-host"
    if fam_p < 0.2 and host_p < 0.4:
        return "polyphyletic / mixed-host"
    return f"{fam_label} / {host_label}"


def summarize_clusters(df: pd.DataFrame, rank: str = "family") -> pd.DataFrame:
    """Summarize clusters at the chosen ICTV rank (one of class/order/family/genus)."""
    if rank not in df.columns:
        raise ValueError(f"rank '{rank}' not in clusters.parquet — "
                         f"available: {[c for c in ('class','order','family','genus') if c in df.columns]}")
    rows = []
    for cl, grp in df.groupby("cluster", observed=True):
        fam_top = _top_label(grp[rank], n=2)
        host_top = _top_label(grp["host_group"], n=2)
        # Exemplars: prefer the dominant rank value (avoids giant-virus bias from
        # raw hit_count), pick 3 most-hit within that group.
        dom_fam = fam_top[0][0] if fam_top else None
        in_dom = grp[grp[rank] == dom_fam] if dom_fam else grp
        if len(in_dom) < 3:
            in_dom = grp
        exemplars = (
            in_dom[["scientific_name", "count"]]
            .dropna(subset=["scientific_name"])
            .sort_values(["count", "scientific_name"], ascending=[False, True])
            .head(3)["scientific_name"]
            .tolist()
        )
        # Flag clusters dominated by a few giant viruses (most hits coming
        # from <10% of cluster members).
        cnts = grp["count"].fillna(0).to_numpy()
        cnts_sorted = sorted(cnts, reverse=True)
        top10pct_n = max(1, len(cnts) // 10)
        top_share = sum(cnts_sorted[:top10pct_n]) / max(sum(cnts), 1)
        giant_flag = "GIANT" if top_share > 0.6 and sum(cnts) > 1000 else ""
        rows.append({
            "cluster": cl,
            "n_viruses": len(grp),
            "n_hits_sum": int(sum(cnts)),
            "top10pct_hit_share": round(top_share, 2),
            "flag": giant_flag,
            "label": _suggest_label(fam_top, host_top),
            f"{rank}_top": _format_share(fam_top),
            "host_top": _format_share(host_top),
            "exemplars": "; ".join(exemplars),
        })
    out = pd.DataFrame(rows)
    # Cluster ids are stringified ints — sort numerically when possible.
    try:
        out["_sort"] = out["cluster"].astype(int)
        out = out.sort_values("n_viruses", ascending=False).drop(columns="_sort")
    except ValueError:
        out = out.sort_values("n_viruses", ascending=False)
    return out


def plot_labeled_umap(df: pd.DataFrame, summary: pd.DataFrame,
                      out_path: Path, footer: str = "",
                      rank: str = "family",
                      marker_summary: pd.DataFrame | None = None,
                      sample_name: str = "") -> None:
    """UMAP scatter colored by cluster, legend shows biological labels.
    If marker_summary is supplied (one row per cluster with 'top_proteins'
    column), appends the top protein names to each legend entry."""
    import matplotlib.pyplot as plt

    # Map cluster id → "c{id}: label" for the legend, ranked by size.
    summary = summary.sort_values("n_viruses", ascending=False).reset_index(drop=True)
    marker_map = {}
    if marker_summary is not None and "top_proteins" in marker_summary.columns:
        marker_map = {
            str(r["cluster"]): str(r.get("top_proteins") or "")
            for _, r in marker_summary.iterrows()
        }
    label_map = {}
    for _, r in summary.iterrows():
        base = f"c{r['cluster']}: {r['label']} (n={r['n_viruses']:,})"
        tp = marker_map.get(str(r["cluster"]), "")
        if tp:
            # Strip [b=...] tags and take top-2 protein names for legend brevity.
            import re as _re
            names = [_re.sub(r"\s*\[b=\d+\]$", "", n.strip())
                     for n in tp.split(";") if n.strip()]
            label_map[str(r["cluster"])] = f"{base}  —  {' / '.join(names[:2])}"
        else:
            label_map[str(r["cluster"])] = base

    TOP_K = 20
    top_ids = [str(c) for c in summary["cluster"].head(TOP_K).tolist()]
    tail_ids = set(str(c) for c in summary["cluster"].iloc[TOP_K:])

    # Decide legend layout: long labels (typical when protein markers attached)
    # go below the plot in 2 columns; short labels stay on the right.
    max_label_len = max((len(v) for v in label_map.values()), default=0)
    legend_below = max_label_len > 55

    if legend_below:
        # Wrap long labels so 2 columns fit cleanly under the plot.
        import textwrap as _tw
        for k, v in list(label_map.items()):
            base, _, proteins = v.partition("  —  ")
            if proteins:
                wrapped_p = "\n    ".join(_tw.wrap(proteins, width=55)) or proteins
                label_map[k] = f"{base}\n    {wrapped_p}"
        fig, ax = plt.subplots(figsize=(14, 10))
    else:
        fig, ax = plt.subplots(figsize=(15, 9))

    # Tail (small clusters) drawn as faint grey backdrop.
    tail_mask = df["cluster"].astype(str).isin(tail_ids)
    if tail_mask.any():
        sub = df[tail_mask]
        n_tail_clust = len(tail_ids)
        ax.scatter(sub["umap_x"], sub["umap_y"], s=6, c="#d8d8d8",
                   alpha=0.4, linewidths=0,
                   label=f"+{n_tail_clust} smaller clusters (n={len(sub):,})")

    # Top-K colored with distinct palette.
    palette = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
    for i, cid in enumerate(top_ids):
        sub = df[df["cluster"].astype(str) == cid]
        if not len(sub):
            continue
        ax.scatter(sub["umap_x"], sub["umap_y"], s=16,
                   c=[palette[i % len(palette)]], alpha=0.85,
                   linewidths=0, label=label_map[cid])
        # Annotate cluster id at centroid for in-plot orientation.
        cx, cy = sub["umap_x"].median(), sub["umap_y"].median()
        ax.annotate(f"c{cid}", (cx, cy), fontsize=9, fontweight="bold",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="black", alpha=0.75, lw=0.5))

    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    title_lead = f"{sample_name}  —  " if sample_name else ""
    ax.set_title(
        f"{title_lead}Leiden clusters with biological labels  —  "
        f"{len(df):,} viruses, {len(summary)} clusters\n"
        f"label = dominant ICTV {rank} / dominant host_group (purity-tagged)"
    )
    if legend_below:
        leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08),
                  fontsize=8, framealpha=0.9, ncol=2,
                  title=f"cluster: {rank} / host  —  top marker proteins",
                  title_fontsize=9, handletextpad=0.5, columnspacing=1.5,
                  borderaxespad=0.3, labelspacing=0.6)
    else:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=8, framealpha=0.9, title=f"cluster: {rank} / host",
                  title_fontsize=9)
    if footer:
        fig.text(0.01, 0.005, footer, fontsize=7, color="gray")
    if legend_below:
        # Don't let tight_layout shrink the axes — pad bottom for the legend instead.
        fig.subplots_adjust(bottom=0.35, left=0.07, right=0.97, top=0.93)
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[label] wrote {out_path}", file=sys.stderr)


def print_markdown(summary: pd.DataFrame, rank: str = "family") -> None:
    print(f"# Cluster labels ({len(summary)} clusters, ICTV rank = {rank})\n")
    print(f"{'cluster':>4}  {'n_v':>5}  {'n_hit':>7}  label  /  {rank}  /  host  /  exemplars")
    print("-" * 110)
    rank_col = f"{rank}_top"
    for _, r in summary.iterrows():
        print(f"{r['cluster']:>4}  {r['n_viruses']:>5}  {r['n_hits_sum']:>7}  "
              f"{r['label']}\n"
              f"           {rank:<6}: {r[rank_col]}\n"
              f"           host:   {r['host_top']}\n"
              f"           ex:     {r['exemplars']}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample-dir", type=Path, required=True,
                    help="Sample output dir containing plots/clusters.parquet")
    ap.add_argument("--rank", default="order",
                    choices=["class", "order", "family", "genus"],
                    help="ICTV rank to label clusters by. Default 'order' matches "
                         "the rank-coherent resolution from scan_resolution.py.")
    ap.add_argument("--with-protein-markers", action="store_true",
                    help="Opt-in: also compute per-cluster marker ORFs and "
                         "annotate them with UniProt protein names. Adds top "
                         "proteins to the labeled UMAP legend and writes "
                         "cluster_markers.tsv + cluster_markers_summary.tsv. "
                         "Default OFF preserves the original analysis output.")
    ap.add_argument("--refs-dir", type=Path,
                    default=Path("/workspace/viral-viz/refs"),
                    help="Used only when --with-protein-markers is set "
                         "(needs refs/cache/uniprot_protein_name.parquet).")
    args = ap.parse_args()

    cl_path = args.sample_dir / "plots" / "clusters.parquet"
    if not cl_path.exists():
        print(f"[label] missing {cl_path} — run plots first", file=sys.stderr)
        sys.exit(1)
    df = pd.read_parquet(cl_path)
    print(f"[label] {len(df):,} viruses across {df['cluster'].nunique()} clusters",
          file=sys.stderr)

    summary = summarize_clusters(df, rank=args.rank)
    out_tsv = args.sample_dir / "plots" / f"cluster_labels_{args.rank}.tsv"
    summary.to_csv(out_tsv, sep="\t", index=False)
    print(f"[label] wrote {out_tsv}", file=sys.stderr)

    marker_summary = None
    if args.with_protein_markers:
        from hvp_viral_viz.cluster_markers import (
            compute_marker_orfs, annotate_markers, summarize as marker_summarize,
            resolve_hits_path, load_hits_for_matrix, load_hits_filtered_to_queries,
        )
        try:
            h_path = resolve_hits_path(args.sample_dir)
        except FileNotFoundError as e:
            print(f"[label]   --with-protein-markers needs hits parquet: {e}",
                  file=sys.stderr)
            sys.exit(1)
        pn_path = args.refs_dir / "cache" / "uniprot_protein_name.parquet"
        pn = pd.read_parquet(pn_path) if pn_path.exists() else None
        if pn is None:
            print(f"[label]   note: {pn_path} missing — protein names will be "
                  f"UniProt-IDs only. Run hvp_viral_viz.uniprot_protein_name first.",
                  file=sys.stderr)
        hits_matrix = load_hits_for_matrix(h_path)
        markers = compute_marker_orfs(hits_matrix, df, n_top=10)
        del hits_matrix
        hits_annot = load_hits_filtered_to_queries(h_path, set(markers["query"]))
        annotated = annotate_markers(markers, hits_annot, pn)
        annotated.to_csv(args.sample_dir / "plots" / "cluster_markers.tsv",
                         sep="\t", index=False)
        marker_summary = marker_summarize(annotated, labels=summary)
        marker_summary.to_csv(args.sample_dir / "plots" / "cluster_markers_summary.tsv",
                              sep="\t", index=False)
        print(f"[label] wrote cluster_markers.tsv + cluster_markers_summary.tsv",
              file=sys.stderr)

    # Need umap coords to plot — these were saved alongside the cluster ids.
    if {"umap_x", "umap_y"}.issubset(df.columns):
        suffix = "_with_proteins" if args.with_protein_markers else ""
        out_png = (args.sample_dir / "plots" /
                   f"umap_by_cluster_labeled_{args.rank}{suffix}.png")
        plot_labeled_umap(df, summary, out_png,
                          footer=f"sample={args.sample_dir.name}  rank={args.rank}",
                          rank=args.rank, marker_summary=marker_summary,
                          sample_name=args.sample_dir.name)

    print_markdown(summary, rank=args.rank)


if __name__ == "__main__":
    main()
