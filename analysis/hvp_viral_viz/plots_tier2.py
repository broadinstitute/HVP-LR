"""Tier 2 visualizations.

Plots:
    1. sunburst_taxonomy      interactive realm→kingdom→phylum→…→family sunburst (plotly HTML)
    2. sankey_orf_to_host     ORF caller → host group flow (matplotlib + plotly)
    3. bits_fident_scatter    bits vs fident hexbin, colored by density; HC overlay
    4. plddt_per_family       boxplots of target avg_plddt for top families
    5. hit_coverage_hist      distribution of (hits per query) and (queries per virus)
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

from hvp_viral_viz.plots import HOST_PALETTE, filter_footer, _save


def plot_sunburst(var: pd.DataFrame, footer: str, out_html: Path) -> None:
    """Interactive plotly sunburst over taxonomy ranks, sized by hit count."""
    import plotly.express as px

    df = var.copy()
    df["count"] = df["count"].astype(int)
    # plotly needs strings, no NaN; AnnData stores str cols as Categorical so coerce first
    for col in ["realm", "kingdom", "phylum", "class", "order", "family"]:
        df[col] = df[col].astype("object").where(df[col].notna(), "unassigned").astype(str)

    fig = px.sunburst(
        df,
        path=["realm", "kingdom", "phylum", "class", "order", "family"],
        values="count",
        title=f"Viral taxonomy sunburst — {df['count'].sum():,} hits across {len(df):,} viruses",
    )
    fig.update_layout(margin=dict(t=60, l=10, r=10, b=40),
                      annotations=[dict(x=0, y=-0.05, xref="paper", yref="paper",
                                        text=footer, showarrow=False, font=dict(size=10, color="gray"))])
    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_html)
    print(f"[plot]   wrote {out_html}", file=sys.stderr)


def plot_sankey_orf_to_host(hits: pd.DataFrame, footer: str, out_html: Path) -> None:
    """Sankey: ORF caller → realm → host_group (hit-count flow)."""
    import plotly.graph_objects as go

    df = hits.copy()
    df["realm"] = df["realm"].astype("object").where(df["realm"].notna(), "unassigned").astype(str)
    df["host_group"] = df["host_group"].astype("object").where(df["host_group"].notna(), "unknown").astype(str)
    df["orf_source"] = df["orf_source"].astype(str)

    # Two flows: orf → realm, realm → host
    f1 = df.groupby(["orf_source", "realm"]).size().reset_index(name="n")
    f2 = df.groupby(["realm", "host_group"]).size().reset_index(name="n")

    nodes = list(dict.fromkeys(
        list(f1["orf_source"]) + list(f1["realm"]) + list(f2["host_group"])
    ))
    idx = {n: i for i, n in enumerate(nodes)}

    src = [idx[s] for s in f1["orf_source"]] + [idx[s] for s in f2["realm"]]
    tgt = [idx[t] for t in f1["realm"]] + [idx[t] for t in f2["host_group"]]
    val = list(f1["n"]) + list(f2["n"])

    # color nodes by category
    node_colors = []
    for n in nodes:
        if n in HOST_PALETTE:
            node_colors.append(HOST_PALETTE[n])
        elif n in {"vs2", "assembly", "genomad", "rescued"}:
            node_colors.append("#444444")
        else:
            node_colors.append("#88aabb")

    fig = go.Figure(go.Sankey(
        node=dict(label=nodes, pad=15, thickness=18, color=node_colors,
                  line=dict(color="black", width=0.3)),
        link=dict(source=src, target=tgt, value=val,
                  color="rgba(150,150,180,0.35)"),
    ))
    fig.update_layout(
        title="ORF caller → realm → host group (foldseek hit counts)",
        margin=dict(t=60, l=10, r=10, b=40),
        annotations=[dict(x=0, y=-0.05, xref="paper", yref="paper",
                          text=footer, showarrow=False, font=dict(size=10, color="gray"))],
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_html)
    print(f"[plot]   wrote {out_html}", file=sys.stderr)


def plot_bits_fident_scatter(hits: pd.DataFrame, footer: str, out: Path) -> None:
    """Hexbin of bits vs fident with high-confidence overlay."""
    fig, ax = plt.subplots(figsize=(8, 7))
    hb = ax.hexbin(
        hits["bits"].astype(float),
        hits["fident"].astype(float),
        gridsize=60,
        cmap="Blues",
        bins="log",
        mincnt=1,
    )
    # Overlay high-conf with red dots (subsample if too many)
    hc = hits[hits["is_high_conf"]]
    if len(hc) > 5000:
        hc = hc.sample(5000, random_state=0)
    ax.scatter(hc["bits"], hc["fident"], s=2, c="red", alpha=0.25, label=f"high-conf (n={int(hits['is_high_conf'].sum()):,})")
    ax.set_xscale("log")
    ax.set_xlabel("bits (log)")
    ax.set_ylabel("fident (sequence identity)")
    ax.set_title("Hit quality: bits × fident (all discovery hits)")
    ax.axvline(300, color="black", lw=0.5, ls="--", alpha=0.6)
    ax.axhline(0.3, color="black", lw=0.5, ls=":", alpha=0.4)
    fig.colorbar(hb, ax=ax, label="hits (log)")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    _save(fig, out, footer)


def plot_plddt_per_family(hits: pd.DataFrame, footer: str, out: Path, top_n: int = 15) -> None:
    """Boxplots of target structure quality (avg_plddt) per family."""
    fam_counts = hits.groupby("family").size().nlargest(top_n).index
    sub = hits[hits["family"].isin(fam_counts)].copy()
    sub["family"] = sub["family"].astype("object")

    # Order by median plddt asc (low quality at top, easier to spot suspect families)
    med = sub.groupby("family")["avg_plddt"].median().sort_values()
    sub["family"] = pd.Categorical(sub["family"], categories=list(med.index), ordered=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * top_n)))
    data = [sub.loc[sub["family"] == f, "avg_plddt"].dropna().values for f in med.index]
    bp = ax.boxplot(data, vert=False, labels=list(med.index), showfliers=False, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#88aacc")
        patch.set_alpha(0.7)
    ax.axvline(70, color="red", lw=0.5, ls="--", alpha=0.6, label="plDDT 70")
    ax.set_xlabel("target avg_plddt (BFVD model quality)")
    ax.set_title(f"Target structure quality per family (top {top_n} by hits)")
    ax.legend(fontsize=8, loc="lower right")
    _save(fig, out, footer)


def plot_hit_coverage_hist(hits: pd.DataFrame, footer: str, out: Path) -> None:
    """Two-panel: (a) hits per query distribution; (b) queries per virus distribution."""
    hits_per_query = hits.groupby("query").size().values
    queries_per_virus = hits.groupby("taxid")["query"].nunique().values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.hist(hits_per_query, bins=np.logspace(0, np.log10(max(hits_per_query) + 1), 40), color="#1f77b4", edgecolor="white")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("hits per query (ORF)")
    ax1.set_ylabel("number of ORFs (log)")
    ax1.set_title(f"Hit fanout per ORF (n={len(hits_per_query):,} ORFs)")
    med = int(np.median(hits_per_query))
    ax1.axvline(med, color="red", lw=0.5, ls="--", label=f"median = {med}")
    ax1.legend(fontsize=8)

    ax2.hist(queries_per_virus, bins=np.logspace(0, np.log10(max(queries_per_virus) + 1), 40),
             color="#2ca02c", edgecolor="white")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("unique ORFs hitting virus")
    ax2.set_ylabel("number of viruses (log)")
    ax2.set_title(f"ORF coverage per virus (n={len(queries_per_virus):,} viruses)")
    med = int(np.median(queries_per_virus))
    ax2.axvline(med, color="red", lw=0.5, ls="--", label=f"median = {med}")
    ax2.legend(fontsize=8)

    _save(fig, out, footer)


def make_all_plots_tier2(sample_dir: Path) -> None:
    plots_dir = sample_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print(f"[plot] loading anndata + hits", file=sys.stderr)
    adata = ad.read_h5ad(sample_dir / "anndata.h5ad")
    var = adata.var.reset_index().rename(columns={"index": "taxid_str"})
    hits = pd.read_parquet(sample_dir / "hits_filtered.parquet")
    footer = filter_footer(adata)

    plot_sunburst(var, footer, plots_dir / "sunburst_taxonomy.html")
    plot_sankey_orf_to_host(hits, footer, plots_dir / "sankey_orf_to_host.html")
    plot_bits_fident_scatter(hits, footer, plots_dir / "bits_fident_scatter.png")
    plot_plddt_per_family(hits, footer, plots_dir / "plddt_per_family.png")
    plot_hit_coverage_hist(hits, footer, plots_dir / "hit_coverage_hist.png")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample-dir", type=Path, required=True)
    args = ap.parse_args()
    make_all_plots_tier2(args.sample_dir)


if __name__ == "__main__":
    main()
