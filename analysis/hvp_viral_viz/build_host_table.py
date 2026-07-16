"""Build per-taxid host annotation table.

Layered fill (per HOST_FILL.md, expanded v2):
  Layer 1 — ICTV VMR Host source, joined hierarchically (species → genus → family).
  Layer 2 — UniProt `virus_hosts` cross-reference (via uniprot_host cache).
            Majority-vote host_group across uniprots-per-taxid.
  Layer 3 — Rule-based by BFVD lineage (class/order/family).
  Layer 4 — Scientific-name regex (phage, NCLDV genera, plant tokens).

Outputs:
    refs/taxid_host.parquet
        taxid, host_group, host_source, host_confidence
    refs/host_fill_audit.tsv  (full table for review)
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

HOST_GROUPS = {
    "human",
    "vertebrate_nonhuman",
    "arthropod",
    "metazoan_other",
    "plant",
    "fungus",
    "protist",
    "bacteria",
    "archaea",
    "unknown",
}

# ---------------------------------------------------------------------------
# Layer 1 — ICTV
# ---------------------------------------------------------------------------

ICTV_VOCAB = {
    "bacteria": "bacteria",
    "archaea": "archaea",
    "fungi": "fungus",
    "plants": "plant",
    "protists": "protist",
    "algae": "protist",
    "vertebrates": "vertebrate_nonhuman",
    "invertebrates": "arthropod",
}


def map_ictv_host(value: str) -> str:
    if not isinstance(value, str):
        return "unknown"
    tokens = [t.strip() for t in re.split(r"[,;]", value) if t.strip()]
    bio_tokens = []
    for t in tokens:
        if "(S)" in t or t.endswith("(S)"):
            continue
        bio_tokens.append(t.lower())
    if not bio_tokens:
        return "unknown"
    for t in bio_tokens:
        if t in ICTV_VOCAB:
            return ICTV_VOCAB[t]
    return "unknown"


def load_ictv(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="VMR MSL41")
    df = df[["Species", "Genus", "Family", "Host source"]].copy()
    df["host_group"] = df["Host source"].map(map_ictv_host)
    return df


def majority_host_per_rank(ictv: pd.DataFrame, rank_col: str) -> pd.DataFrame:
    g = (
        ictv.dropna(subset=[rank_col])
        .groupby([rank_col, "host_group"])
        .size()
        .reset_index(name="n")
    )
    g = g.sort_values([rank_col, "n", "host_group"], ascending=[True, False, True])
    top = g.drop_duplicates(subset=[rank_col], keep="first").rename(columns={rank_col: "name"})
    return top[["name", "host_group", "n"]]


# ---------------------------------------------------------------------------
# Layer 3 — lineage rules
# ---------------------------------------------------------------------------

# Class-level rules — derived from ICTV majorities where frac >= 0.7 OR
# class is monophyletic on the host. Includes rank, name, host_group.
# Generated from ICTV MSL41 (see JOURNAL entry for derivation script).
CLASS_RULES = [
    ("class", "Alsuviricetes", "plant"),
    ("class", "Amabiliviricetes", "fungus"),
    ("class", "Caminiviricetes", "archaea"),
    ("class", "Cardeaviricetes", "vertebrate_nonhuman"),
    ("class", "Caudoviricetes", "bacteria"),
    ("class", "Faserviricetes", "bacteria"),
    ("class", "Herviviricetes", "vertebrate_nonhuman"),
    ("class", "Howeltoviricetes", "fungus"),
    ("class", "Huolimaviricetes", "archaea"),
    ("class", "Laserviricetes", "archaea"),
    ("class", "Leviviricetes", "bacteria"),
    ("class", "Miaviricetes", "fungus"),
    ("class", "Microviricetes", "bacteria"),
    ("class", "Milneviricetes", "plant"),
    ("class", "Mycopleornaviricetes", "fungus"),
    ("class", "Naldaviricetes", "arthropod"),
    ("class", "Papovaviricetes", "vertebrate_nonhuman"),
    ("class", "Pharingeaviricetes", "vertebrate_nonhuman"),
    ("class", "Quintoviricetes", "vertebrate_nonhuman"),
    ("class", "Repensiviricetes", "plant"),
    ("class", "Suforviricetes", "fungus"),
    ("class", "Tectiliviricetes", "bacteria"),
    ("class", "Tokiviricetes", "archaea"),
    ("class", "Tolucaviricetes", "plant"),
    ("class", "Vidaverviricetes", "bacteria"),
]

# Order-level rules — used when class is empty but order is set
ORDER_RULES = [
    ("order", "Ageovirales", "archaea"),
    ("order", "Algavirales", "protist"),
    ("order", "Alpavirales", "bacteria"),
    ("order", "Amoyvirales", "bacteria"),
    ("order", "Asfuvirales", "arthropod"),
    ("order", "Autographivirales", "bacteria"),
    ("order", "Blubervirales", "vertebrate_nonhuman"),
    ("order", "Bullavirales", "bacteria"),
    ("order", "Chitovirales", "vertebrate_nonhuman"),
    ("order", "Cirlivirales", "vertebrate_nonhuman"),
    ("order", "Crassvirales", "bacteria"),
    ("order", "Cryppavirales", "fungus"),
    ("order", "Crytulvirales", "fungus"),
    ("order", "Geplafuvirales", "plant"),
    ("order", "Gokushovirales", "bacteria"),
    ("order", "Grandevirales", "bacteria"),
    ("order", "Haloruvirales", "archaea"),
    ("order", "Herpesvirales", "vertebrate_nonhuman"),
    ("order", "Hypofuvirales", "fungus"),
    ("order", "Imitervirales", "protist"),
    ("order", "Jingchuvirales", "arthropod"),
    ("order", "Kalamavirales", "bacteria"),
    ("order", "Kirjokansivirales", "archaea"),
    ("order", "Lefavirales", "arthropod"),
    ("order", "Ligamenvirales", "archaea"),
    ("order", "Lineavirales", "arthropod"),
    ("order", "Martellivirales", "plant"),
    ("order", "Methanobavirales", "archaea"),
    ("order", "Mindivirales", "bacteria"),
    ("order", "Muvirales", "arthropod"),
    ("order", "Naedrevirales", "plant"),
    ("order", "Nidovirales", "vertebrate_nonhuman"),
    ("order", "Norzivirales", "bacteria"),
    ("order", "Ourlivirales", "fungus"),
    ("order", "Pantevenvirales", "bacteria"),
    ("order", "Patatavirales", "plant"),
    ("order", "Tymovirales", "plant"),
    ("order", "Reovirales", "vertebrate_nonhuman"),
    ("order", "Articulavirales", "vertebrate_nonhuman"),
    ("order", "Bunyavirales", "vertebrate_nonhuman"),  # majority arthropod-borne but most BFVD entries vertebrate
    ("order", "Mononegavirales", "vertebrate_nonhuman"),
    ("order", "Picornavirales", "vertebrate_nonhuman"),
]

# Family-level rules — used when class/order are empty but family is set
FAMILY_RULES = [
    # Bacteria
    ("family", "Microviridae", "bacteria"),
    ("family", "Inoviridae", "bacteria"),
    ("family", "Tectiviridae", "bacteria"),
    ("family", "Plasmaviridae", "bacteria"),
    ("family", "Corticoviridae", "bacteria"),
    ("family", "Sphaerolipoviridae", "bacteria"),
    ("family", "Leviviridae", "bacteria"),
    ("family", "Cystoviridae", "bacteria"),
    ("family", "Autographiviridae", "bacteria"),
    ("family", "Herelleviridae", "bacteria"),
    ("family", "Drexlerviridae", "bacteria"),
    ("family", "Demerecviridae", "bacteria"),
    ("family", "Ackermannviridae", "bacteria"),
    # Archaea
    ("family", "Fuselloviridae", "archaea"),
    ("family", "Rudiviridae", "archaea"),
    ("family", "Lipothrixviridae", "archaea"),
    ("family", "Globuloviridae", "archaea"),
    ("family", "Ampullaviridae", "archaea"),
    ("family", "Bicaudaviridae", "archaea"),
    ("family", "Pleolipoviridae", "archaea"),
    ("family", "Spiraviridae", "archaea"),
    ("family", "Salterprovirus", "archaea"),
    # Plant
    ("family", "Geminiviridae", "plant"),
    ("family", "Caulimoviridae", "plant"),
    ("family", "Bromoviridae", "plant"),
    ("family", "Closteroviridae", "plant"),
    ("family", "Tymoviridae", "plant"),
    ("family", "Tobamoviridae", "plant"),
    ("family", "Virgaviridae", "plant"),
    ("family", "Tombusviridae", "plant"),
    ("family", "Luteoviridae", "plant"),
    ("family", "Potyviridae", "plant"),
    ("family", "Secoviridae", "plant"),
    # Protist (giant viruses, algal)
    ("family", "Mimiviridae", "protist"),
    ("family", "Marseilleviridae", "protist"),
    ("family", "Phycodnaviridae", "protist"),
    ("family", "Pithoviridae", "protist"),
    ("family", "Pandoraviridae", "protist"),
    # Vertebrate
    ("family", "Iridoviridae", "vertebrate_nonhuman"),
    ("family", "Asfarviridae", "vertebrate_nonhuman"),
    ("family", "Polyomaviridae", "vertebrate_nonhuman"),
    ("family", "Papillomaviridae", "vertebrate_nonhuman"),
    ("family", "Herpesviridae", "vertebrate_nonhuman"),
    ("family", "Coronaviridae", "vertebrate_nonhuman"),
    ("family", "Filoviridae", "vertebrate_nonhuman"),
    ("family", "Retroviridae", "vertebrate_nonhuman"),
    ("family", "Flaviviridae", "vertebrate_nonhuman"),
    ("family", "Paramyxoviridae", "vertebrate_nonhuman"),
    ("family", "Orthomyxoviridae", "vertebrate_nonhuman"),
    ("family", "Adenoviridae", "vertebrate_nonhuman"),
    ("family", "Reoviridae", "vertebrate_nonhuman"),
    ("family", "Picornaviridae", "vertebrate_nonhuman"),
    ("family", "Hepadnaviridae", "vertebrate_nonhuman"),
    ("family", "Poxviridae", "vertebrate_nonhuman"),
    ("family", "Rhabdoviridae", "vertebrate_nonhuman"),
    ("family", "Bunyaviridae", "vertebrate_nonhuman"),
    ("family", "Arenaviridae", "vertebrate_nonhuman"),
    ("family", "Anelloviridae", "vertebrate_nonhuman"),
    ("family", "Astroviridae", "vertebrate_nonhuman"),
    ("family", "Caliciviridae", "vertebrate_nonhuman"),
    ("family", "Togaviridae", "vertebrate_nonhuman"),
    ("family", "Hepeviridae", "vertebrate_nonhuman"),
    ("family", "Parvoviridae", "vertebrate_nonhuman"),
    ("family", "Circoviridae", "vertebrate_nonhuman"),
    # Arthropod
    ("family", "Baculoviridae", "arthropod"),
    ("family", "Polydnaviridae", "arthropod"),
    ("family", "Nudiviridae", "arthropod"),
    ("family", "Hytrosaviridae", "arthropod"),
    ("family", "Ascoviridae", "arthropod"),
    ("family", "Dicistroviridae", "arthropod"),
    # Fungal
    ("family", "Hypoviridae", "fungus"),
    ("family", "Partitiviridae", "fungus"),
    ("family", "Totiviridae", "fungus"),
    ("family", "Narnaviridae", "fungus"),
]

# Layer 3 lineage rules — order matters: family first (most specific),
# then order, then class. This keeps prior specificity-first behavior.
LINEAGE_RULES = FAMILY_RULES + ORDER_RULES + CLASS_RULES


# ---------------------------------------------------------------------------
# Layer 4 — name regex
# ---------------------------------------------------------------------------

PHAGE_RE = re.compile(r"\b(phage|prophage|bacteriophage)\b", re.IGNORECASE)
PROK_DNA_RE = re.compile(r"\bProkaryotic\b.+\bvirus\b", re.IGNORECASE)

# NCLDV / giant virus genera (Acanthamoeba-infecting, etc.)
NCLDV_GENERA = (
    r"Pithovirus|Pandoravirus|Pacmanvirus|Faustovirus|Kaumoebavirus|"
    r"Tupanvirus|Mollivirus|Cedratvirus|Orpheovirus|Yaravirus|Medusavirus|"
    r"Klosneuvirus|Catovirus|Hokovirus|Indivirus|Mimivirus|Megavirus|"
    r"Marseillevirus|Melbournevirus|Lausannevirus|Cannes 8 virus|"
    r"Tokyovirus|Clandestinovirus|Sylvanvirus|Solumvirus"
)
NCLDV_RE = re.compile(r"\b(?:" + NCLDV_GENERA + r")\b")

PLANT_RE = re.compile(
    r"\b(tobacco|tomato|potato|maize|wheat|barley|rice|cassava|cucumber|pepper|grape|"
    r"banana|citrus|cotton|lettuce|onion|pea|bean|soybean|squash|sugarcane|"
    r"watermelon|carnation|chrysanthemum|tulip|orchid|rose|"
    r"clover|trifolium|bromus|ash|mosaic)\b",
    re.IGNORECASE,
)


def name_regex_host(name: str) -> str | None:
    if not isinstance(name, str) or not name:
        return None
    if PROK_DNA_RE.search(name):
        return "bacteria"
    if PHAGE_RE.search(name):
        return "bacteria"
    if NCLDV_RE.search(name):
        return "protist"
    if PLANT_RE.search(name):
        return "plant"
    return None


# ---------------------------------------------------------------------------
# Human override taxid list
# ---------------------------------------------------------------------------

HUMAN_TAXIDS = {
    2697049, 694009, 1335626,
    11137, 31631, 277944, 290028,
    11676, 11709,
    10298, 10310,
    10335,
    10376,
    10359,
    32603, 32604, 10372, 37296,
    10407, 11103, 12475, 12092, 1678143,
    333760, 333761, 10566, 333922, 337051, 337041,
    11320, 11520, 11552,
    11250, 12814, 208893, 208895,
    11234, 11161, 11041,
    11983, 95341, 1216201,
    28875,
    10535, 28285, 28286, 130309, 565302,
    186538, 565995,
    11269, 33727,
    11053, 11060, 11069, 11070,
    64320, 11082, 11089,
    12080, 12083, 12086,
    11292,
    10632, 1891762, 10798,
    10255, 10245, 10244, 10243,
}


# ---------------------------------------------------------------------------
# Layer 2 — UniProt virus_hosts → per-taxid majority
# ---------------------------------------------------------------------------

def uniprot_host_per_taxid(
    uniprot_taxid: pd.DataFrame,
    uniprot_host: pd.DataFrame,
    host_taxid_to_group: pd.DataFrame,
) -> pd.DataFrame:
    """For each viral taxid, derive majority host_group from its uniprots.

    Pipeline per viral taxid:
        uniprots-of-taxid → host_taxid list (multi-host entries split) →
        host_group via taxid_to_group → majority vote across all uniprots.

    Returns: DataFrame with columns [taxid, host_group_uniprot, n_uniprots_with_host]
    """
    # Build host_taxid → host_group lookup
    h2g = dict(zip(host_taxid_to_group["host_taxid"].astype(str), host_taxid_to_group["host_group"]))

    # Expand uniprot_host's host_taxid (';'-joined) into one row per (uniprot, host_taxid)
    uh = uniprot_host[uniprot_host["host_taxid"].astype(str).str.len() > 0].copy()
    uh["host_taxid"] = uh["host_taxid"].astype(str).str.split(";")
    uh = uh.explode("host_taxid")
    uh = uh[uh["host_taxid"].str.len() > 0]
    uh["host_group"] = uh["host_taxid"].map(h2g).fillna("unknown")
    uh = uh[uh["host_group"] != "unknown"]

    # Join uniprot → viral taxid (uniprot_taxid table)
    joined = uniprot_taxid.merge(uh[["uniprot_id", "host_group"]], on="uniprot_id", how="inner")

    # Majority host_group per viral taxid (ties → alphabetic for determinism)
    counts = joined.groupby(["taxid", "host_group"]).size().reset_index(name="n")
    counts = counts.sort_values(["taxid", "n", "host_group"], ascending=[True, False, True])
    top = counts.drop_duplicates(subset=["taxid"], keep="first")

    # Coverage: number of distinct uniprots with any host annotation per taxid
    nuni = joined.groupby("taxid")["uniprot_id"].nunique().reset_index(name="n_uniprots_with_host")

    out = top.merge(nuni, on="taxid").rename(columns={"host_group": "host_group_uniprot"})
    return out[["taxid", "host_group_uniprot", "n_uniprots_with_host"]]


# ---------------------------------------------------------------------------
# Main build
# ---------------------------------------------------------------------------

def build_host_table(refs_dir: Path, ictv_path: Path) -> None:
    print(f"[host] loading BFVD lineage", file=sys.stderr)
    lineage = pd.read_parquet(refs_dir / "taxid_lineage.parquet")
    print(f"[host]   {len(lineage):,} taxids", file=sys.stderr)

    print(f"[host] loading ICTV VMR from {ictv_path}", file=sys.stderr)
    ictv = load_ictv(ictv_path)
    print(f"[host]   {len(ictv):,} ICTV rows", file=sys.stderr)

    species_lookup = majority_host_per_rank(ictv, "Species")
    genus_lookup = majority_host_per_rank(ictv, "Genus")
    family_lookup = majority_host_per_rank(ictv, "Family")
    print(
        f"[host]   ICTV lookups: species={len(species_lookup):,}, "
        f"genus={len(genus_lookup):,}, family={len(family_lookup):,}",
        file=sys.stderr,
    )

    # ---- Layer 1: hierarchical ICTV ----
    out = lineage.copy()
    out["host_group"] = pd.NA
    out["host_source"] = pd.NA

    sp = species_lookup.set_index("name")["host_group"]
    mask = out["species"].isin(sp.index) & out["host_group"].isna()
    out.loc[mask, "host_group"] = out.loc[mask, "species"].map(sp)
    out.loc[mask, "host_source"] = "ictv_species"
    n_sp = mask.sum()

    gn = genus_lookup.set_index("name")["host_group"]
    mask = out["genus"].isin(gn.index) & out["host_group"].isna()
    out.loc[mask, "host_group"] = out.loc[mask, "genus"].map(gn)
    out.loc[mask, "host_source"] = "ictv_genus"
    n_gn = mask.sum()

    fm = family_lookup.set_index("name")["host_group"]
    mask = out["family"].isin(fm.index) & out["host_group"].isna()
    out.loc[mask, "host_group"] = out.loc[mask, "family"].map(fm)
    out.loc[mask, "host_source"] = "ictv_family"
    n_fm = mask.sum()

    print(
        f"[host] L1 ICTV: species={n_sp:,}, genus={n_gn:,}, family={n_fm:,} "
        f"(cumulative {n_sp + n_gn + n_fm:,}/{len(out):,})",
        file=sys.stderr,
    )

    # ---- Layer 2: UniProt virus_hosts ----
    uh_cache = refs_dir / "cache" / "uniprot_host.parquet"
    h2g_cache = refs_dir / "cache" / "host_taxid_to_group.parquet"
    n_uniprot = 0
    if uh_cache.exists() and h2g_cache.exists():
        uh = pd.read_parquet(uh_cache)
        h2g = pd.read_parquet(h2g_cache)
        ut = pd.read_parquet(refs_dir / "uniprot_taxid.parquet")
        uphost = uniprot_host_per_taxid(ut, uh, h2g)
        # Apply only to rows still unfilled
        unfilled_mask = out["host_group"].isna()
        merged = out.loc[unfilled_mask, ["taxid"]].merge(uphost, on="taxid", how="inner")
        fill_idx = out.index[out["taxid"].isin(set(merged["taxid"]))]
        # Build map taxid -> host_group_uniprot for unfilled set
        map_uni = dict(zip(merged["taxid"], merged["host_group_uniprot"]))
        for i in fill_idx:
            if pd.isna(out.at[i, "host_group"]):
                hg = map_uni.get(out.at[i, "taxid"])
                if hg:
                    out.at[i, "host_group"] = hg
                    out.at[i, "host_source"] = "uniprot_host"
                    n_uniprot += 1
        print(f"[host] L2 UniProt: filled {n_uniprot:,}", file=sys.stderr)
    else:
        print(f"[host] L2 UniProt: caches missing (run uniprot_host + host_taxid_lookup); skipping", file=sys.stderr)

    # ---- Layer 3: lineage rules ----
    n_rule = 0
    for rank_col, name, hg in LINEAGE_RULES:
        mask = (out[rank_col] == name) & out["host_group"].isna()
        n = mask.sum()
        if n:
            out.loc[mask, "host_group"] = hg
            out.loc[mask, "host_source"] = f"rule_{rank_col}"
            n_rule += n
    print(f"[host] L3 rules: filled {n_rule:,}", file=sys.stderr)

    # ---- Layer 4: name regex ----
    mask = out["host_group"].isna()
    name_fill = out.loc[mask, "scientific_name"].map(name_regex_host)
    out.loc[mask, "host_group"] = out.loc[mask, "host_group"].fillna(name_fill)
    n_name_filled = mask.sum() - out["host_group"].isna().sum()
    out.loc[mask & out["host_source"].isna() & out["host_group"].notna(), "host_source"] = "name_regex"
    print(f"[host] L4 name regex: filled {n_name_filled:,}", file=sys.stderr)

    # ---- Final: unknowns ----
    mask = out["host_group"].isna()
    out.loc[mask, "host_group"] = "unknown"
    out.loc[mask, "host_source"] = "unfilled"
    print(f"[host] unfilled → unknown: {mask.sum():,}", file=sys.stderr)

    # ---- Human override ----
    mask = out["taxid"].isin(HUMAN_TAXIDS)
    n_human = mask.sum()
    out.loc[mask, "host_group"] = "human"
    out.loc[mask, "host_source"] = out.loc[mask, "host_source"].astype(str) + "+human_override"
    print(f"[host] human override: {n_human:,} taxids", file=sys.stderr)

    # ---- host_confidence ----
    conf_map = {
        "ictv_species": "high",
        "ictv_genus": "high",
        "ictv_family": "medium",
        "uniprot_host": "high",
        "rule_class": "medium",
        "rule_order": "medium",
        "rule_family": "medium",
        "name_regex": "low",
        "unfilled": "low",
    }

    def conf(src):
        if not isinstance(src, str):
            return "low"
        base = src.split("+", 1)[0]
        return conf_map.get(base, "low")

    out["host_confidence"] = out["host_source"].map(conf)

    # ---- Summary ----
    print("\n[host] host_group distribution:", file=sys.stderr)
    print(out["host_group"].value_counts().to_string(), file=sys.stderr)
    print("\n[host] host_source distribution:", file=sys.stderr)
    print(out["host_source"].value_counts().to_string(), file=sys.stderr)

    main_cols = ["taxid", "host_group", "host_source", "host_confidence"]
    audit_cols = ["taxid", "scientific_name", "family", "genus", "species",
                  "host_group", "host_source", "host_confidence"]

    out[main_cols].to_parquet(refs_dir / "taxid_host.parquet", index=False)
    print(f"[host]   wrote {refs_dir / 'taxid_host.parquet'}", file=sys.stderr)
    out[audit_cols].to_csv(refs_dir / "host_fill_audit.tsv", sep="\t", index=False)
    print(f"[host]   wrote {refs_dir / 'host_fill_audit.tsv'}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument(
        "--ictv-vmr",
        type=Path,
        default=Path("/workspace/viral-viz/refs/downloads/ictv_vmr_msl41.xlsx"),
    )
    args = ap.parse_args()
    build_host_table(args.refs_dir, args.ictv_vmr)


if __name__ == "__main__":
    main()
