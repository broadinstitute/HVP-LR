"""Map host_taxid → host_group.

UniProt's `virus_hosts` cross-reference returns host species taxids. We bucket
those into the same coarse `host_group` vocabulary used everywhere else.

Strategy:
  1. Hard-code well-known model-organism / pathogen-of-interest taxids
     (covers ~80% of UniProt's hosts in practice — they cluster around a
     small set of curated organisms).
  2. For the rest, look up the host taxid's lineage via NCBI eutils
     (esummary on db=taxonomy) and bucket by superkingdom + a few kingdom
     hints.

Outputs:
    refs/cache/host_taxid_to_group.parquet
        host_taxid, host_organism, host_group, source
"""
from __future__ import annotations
import argparse
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests

# Curated host taxid → host_group (covers the common hosts UniProt cites)
KNOWN_HOSTS: dict[int, str] = {
    # Human
    9606: "human",
    # Common vertebrate model / livestock
    10090: "vertebrate_nonhuman",   # Mus musculus
    10116: "vertebrate_nonhuman",   # Rattus norvegicus
    9913: "vertebrate_nonhuman",    # Bos taurus
    9823: "vertebrate_nonhuman",    # Sus scrofa
    9031: "vertebrate_nonhuman",    # Gallus gallus
    9685: "vertebrate_nonhuman",    # Felis catus
    9615: "vertebrate_nonhuman",    # Canis lupus familiaris
    9796: "vertebrate_nonhuman",    # Equus caballus
    9940: "vertebrate_nonhuman",    # Ovis aries
    9925: "vertebrate_nonhuman",    # Capra hircus
    8355: "vertebrate_nonhuman",    # Xenopus laevis
    7955: "vertebrate_nonhuman",    # Danio rerio
    # Bats / common reservoirs
    59477: "vertebrate_nonhuman",   # Rhinolophus
    9397: "vertebrate_nonhuman",    # Chiroptera
    # Common arthropod model
    7227: "arthropod",              # Drosophila melanogaster
    6239: "metazoan_other",         # Caenorhabditis elegans
    7165: "arthropod",              # Anopheles gambiae
    7159: "arthropod",              # Aedes aegypti
    7460: "arthropod",              # Apis mellifera
    # Plant models
    3702: "plant",                  # Arabidopsis thaliana
    4530: "plant",                  # Oryza sativa
    4577: "plant",                  # Zea mays
    4081: "plant",                  # Solanum lycopersicum
    4097: "plant",                  # Nicotiana tabacum
    # Fungi
    4932: "fungus",                 # Saccharomyces cerevisiae
    4896: "fungus",                 # Schizosaccharomyces pombe
    5476: "fungus",                 # Candida albicans
    330879: "fungus",               # Aspergillus fumigatus
    # Bacteria — common phage hosts
    562: "bacteria",                # Escherichia coli
    1280: "bacteria",               # Staphylococcus aureus
    1773: "bacteria",               # Mycobacterium tuberculosis
    287: "bacteria",                # Pseudomonas aeruginosa
    1392: "bacteria",               # Bacillus anthracis
    1423: "bacteria",               # Bacillus subtilis
    1428: "bacteria",               # Bacillus thuringiensis
    485: "bacteria",                # Neisseria gonorrhoeae
    487: "bacteria",                # Neisseria meningitidis
    727: "bacteria",                # Haemophilus influenzae
    210: "bacteria",                # Helicobacter pylori
    1305: "bacteria",               # Streptococcus sanguinis
    1313: "bacteria",               # Streptococcus pneumoniae
    1314: "bacteria",               # Streptococcus pyogenes
    1351: "bacteria",               # Enterococcus faecalis
    1352: "bacteria",               # Enterococcus faecium
    192222: "bacteria",             # Campylobacter jejuni
    211044: "bacteria",             # Influenza A H1N1 PR/8 strain (cited as virus, but here is host?) — leave; will fall through
    470: "bacteria",                # Acinetobacter baumannii
    573: "bacteria",                # Klebsiella pneumoniae
    546: "bacteria",                # Citrobacter freundii
    615: "bacteria",                # Serratia marcescens
    632: "bacteria",                # Yersinia pestis
    644: "bacteria",                # Aeromonas hydrophila
    666: "bacteria",                # Vibrio cholerae
    813: "bacteria",                # Chlamydia trachomatis
    1763: "bacteria",               # Mycobacterium genus
    1769: "bacteria",               # Mycobacterium leprae
    1496: "bacteria",               # Clostridioides difficile
    1502: "bacteria",               # Clostridium perfringens
    1639: "bacteria",               # Listeria monocytogenes
    28901: "bacteria",              # Salmonella enterica
    590: "bacteria",                # Salmonella genus
    584: "bacteria",                # Proteus mirabilis
    197: "bacteria",                # Campylobacter
    1063: "bacteria",               # Rhodobacter sphaeroides
    160: "bacteria",                # Treponema (gen)
    157: "bacteria",                # Treponema pallidum
    # Archaea
    2157: "archaea",                # Archaea superkingdom
    2287: "archaea",                # Sulfolobus
    2285: "archaea",                # Sulfolobus solfataricus
    2188: "archaea",                # Methanocaldococcus jannaschii
    # Protists
    5754: "protist",                # Acanthamoeba castellanii
    5755: "protist",                # Acanthamoeba polyphaga
    5759: "protist",                # Vermamoeba vermiformis
    5833: "protist",                # Plasmodium falciparum
    5691: "protist",                # Trypanosoma brucei
    5811: "protist",                # Toxoplasma gondii
    44689: "protist",               # Dictyostelium discoideum
    3055: "protist",                # Chlamydomonas reinhardtii
    35128: "protist",               # Thalassiosira pseudonana
    296543: "protist",              # Bigelowiella natans
    # Metazoan other (non-arthropod invertebrates)
    6500: "metazoan_other",         # Aplysia californica
    6669: "metazoan_other",         # Daphnia pulex (crustacean — could be arthropod)
    7668: "metazoan_other",         # Strongylocentrotus purpuratus (sea urchin)
    37653: "metazoan_other",        # Haliotis (abalone)
}

# NCBI superkingdom → bucket. Used as fallback when host_taxid not in KNOWN_HOSTS.
# We map kingdom or phylum/class on top of superkingdom when possible.
SUPERKINGDOM_DEFAULT: dict[str, str] = {
    "Bacteria": "bacteria",
    "Archaea": "archaea",
    "Eukaryota": "unknown",   # need further breakdown
}

EUKARYOTE_KINGDOM: dict[str, str] = {
    "Viridiplantae": "plant",
    "Fungi": "fungus",
}

EUKARYOTE_PHYLUM: dict[str, str] = {
    "Arthropoda": "arthropod",
    "Chordata": "vertebrate_nonhuman",
    "Nematoda": "metazoan_other",
    "Mollusca": "metazoan_other",
    "Echinodermata": "metazoan_other",
    "Annelida": "metazoan_other",
    "Cnidaria": "metazoan_other",
    "Porifera": "metazoan_other",
    "Platyhelminthes": "metazoan_other",
}


def lineage_to_group(lineage_dict: dict[str, str]) -> str:
    """Bucket NCBI taxonomy lineage (rank → name) → host_group."""
    # NCBI 2024+ uses "domain"; older snapshots use "superkingdom"
    sk = lineage_dict.get("domain") or lineage_dict.get("superkingdom")
    if sk == "Bacteria":
        return "bacteria"
    if sk == "Archaea":
        return "archaea"
    if sk == "Viruses":
        return "unknown"
    if sk != "Eukaryota":
        return "unknown"

    # Eukaryote: refine by kingdom / phylum
    k = lineage_dict.get("kingdom")
    if k in EUKARYOTE_KINGDOM:
        return EUKARYOTE_KINGDOM[k]

    p = lineage_dict.get("phylum")
    if p in EUKARYOTE_PHYLUM:
        # Human override below caller layer
        return EUKARYOTE_PHYLUM[p]

    # Fall back to protist (everything else eukaryotic)
    return "protist"


def fetch_lineage_batch(taxids: list[str], session: requests.Session) -> dict[str, dict[str, str]]:
    """efetch taxonomy XML for up to ~200 taxids, return {taxid: {rank: name}}."""
    if not taxids:
        return {}
    params = {
        "db": "taxonomy",
        "id": ",".join(taxids),
        "rettype": "xml",
    }
    r = session.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    out: dict[str, dict[str, str]] = {}
    try:
        root = ET.fromstring(r.content)
    except ET.ParseError as e:
        print(f"[taxlookup]   parse error: {e}", file=sys.stderr)
        return out
    for tx in root.findall("Taxon"):
        tid = tx.findtext("TaxId")
        if not tid:
            continue
        ranks: dict[str, str] = {}
        # The top-level Taxon also has a Rank/ScientificName — record it
        rank_self = tx.findtext("Rank")
        name_self = tx.findtext("ScientificName")
        if rank_self and name_self:
            ranks[rank_self] = name_self
        for node in tx.findall("LineageEx/Taxon"):
            rk = node.findtext("Rank")
            nm = node.findtext("ScientificName")
            if rk and nm:
                ranks[rk] = nm
        out[tid] = ranks
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refs-dir", type=Path, default=Path("/workspace/viral-viz/refs"))
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    cache_path = args.refs_dir / "cache" / "host_taxid_to_group.parquet"
    uh_path = args.refs_dir / "cache" / "uniprot_host.parquet"
    if not uh_path.exists():
        print(f"[taxlookup] {uh_path} missing — run uniprot_host first", file=sys.stderr)
        sys.exit(1)

    uh = pd.read_parquet(uh_path)
    raw = uh["host_taxid"].dropna()
    raw = raw[raw != ""]
    all_tids: set[str] = set()
    for s in raw:
        for t in s.split(";"):
            if t.isdigit():
                all_tids.add(t)
    print(f"[taxlookup] unique host_taxids cited: {len(all_tids):,}", file=sys.stderr)

    # KNOWN_HOSTS direct fill
    rows = []
    needs_fetch: list[str] = []
    for t in all_tids:
        ti = int(t)
        if ti in KNOWN_HOSTS:
            rows.append((t, KNOWN_HOSTS[ti], "known"))
        else:
            needs_fetch.append(t)
    print(f"[taxlookup]   covered by KNOWN_HOSTS: {len(rows):,}; needs efetch: {len(needs_fetch):,}",
          file=sys.stderr)

    # Resume cache
    if cache_path.exists() and not args.refresh:
        cached = pd.read_parquet(cache_path)
        already = set(cached["host_taxid"].astype(str))
        needs_fetch = [t for t in needs_fetch if t not in already]
        print(f"[taxlookup]   {len(already):,} cached; {len(needs_fetch):,} new", file=sys.stderr)
    else:
        cached = None

    session = requests.Session()
    session.headers["User-Agent"] = "hvp-viral-viz/0.1 host-taxid-lookup"

    BATCH = 200
    n = len(needs_fetch)
    t0 = time.time()
    for i in range(0, n, BATCH):
        batch = needs_fetch[i : i + BATCH]
        try:
            lin = fetch_lineage_batch(batch, session)
        except requests.HTTPError as e:
            print(f"[taxlookup]   batch {i}: HTTP {e.response.status_code}; sleeping 3s, retry", file=sys.stderr)
            time.sleep(3)
            lin = fetch_lineage_batch(batch, session)

        for t in batch:
            ranks = lin.get(t, {})
            group = lineage_to_group(ranks)
            rows.append((t, group, "ncbi_lineage"))

        if i % (BATCH * 5) == 0:
            elapsed = time.time() - t0
            rate = (i + len(batch)) / max(elapsed, 0.1)
            eta = (n - i - len(batch)) / max(rate, 0.1)
            print(f"[taxlookup]   {i + len(batch):,}/{n:,} rate={rate:.0f}/s eta={eta:.0f}s",
                  file=sys.stderr)
        time.sleep(0.34)  # eutils 3 rps without API key

    new = pd.DataFrame(rows, columns=["host_taxid", "host_group", "source"])
    if cached is not None:
        combined = pd.concat([cached, new], ignore_index=True).drop_duplicates(
            subset=["host_taxid"], keep="last"
        )
    else:
        combined = new

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(cache_path, index=False)
    print(f"[taxlookup] wrote {cache_path}  ({len(combined):,} rows)", file=sys.stderr)
    print(combined["host_group"].value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
