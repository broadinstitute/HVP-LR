version 1.0

import "../../structs/Structs.wdl"

# ── BinSummary ────────────────────────────────────────────────────────────────

task BinSummary {

    meta {
        description: "Produce a one-row-per-bin summary combining skani taxonomy and CheckM2 quality. Selects the best GTDB reference per bin by skani_score = ANI * Align_fraction_query / 100. Computes MIMAG quality tier (high-quality: completeness >= 90 and contamination < 5; medium-quality: completeness >= 50 and contamination < 10; low-quality: otherwise; contaminated: contamination >= 10) and flags single-contig bins."
        tool:         "Python 3"
        outputs: {
            bin_summary_tsv: "Per-bin TSV: bin_name, skani_score, ANI, best_ref_accession, gtdb_taxonomy, completeness, contamination, quality_score, mimag_quality, single_contig, genome_size, contig_n50, total_contigs"
        }
    }

    parameter_meta {
        checkm2_quality_tsv:   "CheckM2 quality_report.tsv from CheckM2 task"
        skani_annotated_tsv:   "Skani annotated results TSV from SkaniAnnotate task"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   checkm2_quality_tsv
        File   skani_annotated_tsv
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python3 - <<'PYEOF'
import csv, os, sys
from collections import Counter

skani_file   = "~{skani_annotated_tsv}"
checkm2_file = "~{checkm2_quality_tsv}"
out_file     = "~{sample_name}.bin_summary.tsv"


def bin_name_from_path(path):
    base = os.path.basename(path)
    if base.endswith(".fa"):
        base = base[:-3]
    return base


def mimag_quality(comp, cont):
    if cont >= 10:
        return "contaminated"
    if comp >= 90 and cont < 5:
        return "high-quality"
    if comp >= 50 and cont < 10:
        return "medium-quality"
    return "low-quality"


def accession_from_ref(ref_path):
    base = os.path.basename(ref_path)
    if base.endswith("_genomic.fna.gz"):
        base = base[: -len("_genomic.fna.gz")]
    elif base.endswith(".gz"):
        base = base[:-3]
    return base


# Pass 1: skani — best hit per bin
best = {}
with open(skani_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row.get("Ref_file") == "Ref_file":
            continue
        name = bin_name_from_path(row["Query_file"])
        try:
            score = float(row["ANI"]) * float(row["Align_fraction_query"]) / 100.0
        except ValueError:
            continue
        if name not in best or score > best[name]["_score"]:
            row["_score"]     = score
            row["_accession"] = accession_from_ref(row["Ref_file"])
            best[name] = row

# Pass 2: CheckM2
checkm2 = {}
with open(checkm2_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        checkm2[row["Name"]] = row

out_cols = [
    "bin_name", "skani_score", "ANI", "Align_fraction_ref", "Align_fraction_query",
    "best_ref_accession", "gtdb_taxonomy",
    "completeness", "contamination", "quality_score", "mimag_quality",
    "single_contig", "genome_size", "contig_n50", "total_contigs",
]

with open(out_file, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(out_cols)
    for name in sorted(checkm2.keys()):
        cm = checkm2[name]
        try:
            comp   = float(cm["Completeness"])
            cont   = float(cm["Contamination"])
            nctg   = int(cm["Total_Contigs"])
            qscore = round(comp - 5.0 * cont, 2)
            tier   = mimag_quality(comp, cont)
            single = 1 if nctg == 1 else 0
        except (ValueError, KeyError):
            comp = cont = nctg = qscore = "NA"
            tier = single = "NA"

        if name in best:
            sk = best[name]
            w.writerow([
                name, round(sk["_score"], 2), sk["ANI"],
                sk["Align_fraction_ref"], sk["Align_fraction_query"],
                sk["_accession"], sk["gtdb_taxonomy"],
                comp, cont, qscore, tier, single,
                cm.get("Genome_Size", "NA"), cm.get("Contig_N50", "NA"), cm.get("Total_Contigs", "NA"),
            ])
        else:
            w.writerow([
                name, "NA", "NA", "NA", "NA", "NA", "NA",
                comp, cont, qscore, tier, single,
                cm.get("Genome_Size", "NA"), cm.get("Contig_N50", "NA"), cm.get("Total_Contigs", "NA"),
            ])

n_bins  = len(checkm2)
n_skani = len(best)
print(f"Wrote {n_bins} bins ({n_skani} with skani hit, {n_bins - n_skani} with no hit) -> {out_file}", file=sys.stderr)

tier_counts = Counter()
with open(out_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        tier_counts[row["mimag_quality"]] += 1
for tier in ["high-quality", "medium-quality", "low-quality", "contaminated", "NA"]:
    if tier_counts[tier]:
        print(f"  {tier}: {tier_counts[tier]}", file=sys.stderr)
PYEOF
    >>>

    output {
        File bin_summary_tsv = "~{sample_name}.bin_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            4,
        disk_gb:           20,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:          select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:       select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb,  default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:   select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:      select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:       select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:           select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# ── MagSummary ────────────────────────────────────────────────────────────────

task MagSummary {

    meta {
        description: "Aggregate bin-level metrics from BinSummary into per-sample scalar outputs for the Terra data table. Also writes a combined TSV joining assembly stats with bin aggregates. Scalar outputs follow the same read_int()/read_float() pattern as ViralOverallSummary."
        tool:         "Python 3"
        outputs: {
            mag_stats_tsv:             "Single-row TSV of all assembly and bin aggregate metrics",
            total_bins:                "Total number of DAS_Tool-refined bins",
            n_high_quality:            "Bins with completeness >= 90 and contamination < 5 (MIMAG HQ)",
            n_medium_quality:          "Bins with completeness >= 50 and contamination < 10 (MIMAG MQ)",
            n_low_quality:             "Bins with completeness < 50 and contamination < 10 (MIMAG LQ)",
            n_contaminated:            "Bins with contamination >= 10 (MIMAG contaminated)",
            n_no_skani_hit:            "Bins with no skani reference hit",
            n_single_contig:           "Bins composed of a single contig (potential complete circular MAG)",
            n_single_contig_hq:        "HQ bins composed of a single contig",
            total_bases_all_bins:      "Total base pairs across all bins",
            total_bases_hq:            "Total base pairs in HQ bins",
            total_bases_medium:        "Total base pairs in medium-quality bins",
            pct_assembly_bases_binned: "Percentage of assembly bases captured in any bin",
            pct_assembly_bases_hq:     "Percentage of assembly bases captured in HQ bins",
            mean_completeness_all:     "Mean CheckM2 completeness across all bins",
            mean_completeness_hq:      "Mean CheckM2 completeness across HQ bins",
            mean_contamination_all:    "Mean CheckM2 contamination across all bins",
            mean_contamination_hq:     "Mean CheckM2 contamination across HQ bins",
            mean_quality_score_hq:     "Mean CheckM2 quality score (completeness - 5*contamination) across HQ bins",
            n_distinct_species_all:    "Number of distinct GTDB species-level assignments across all bins",
            n_distinct_species_hq:     "Number of distinct GTDB species-level assignments among HQ bins",
            n_distinct_genera_all:     "Number of distinct GTDB genus-level assignments across all bins",
            n_distinct_genera_hq:      "Number of distinct GTDB genus-level assignments among HQ bins"
        }
    }

    parameter_meta {
        bin_summary_tsv:       "Per-bin summary TSV from BinSummary"
        asm_stats_tsv:         "Assembly seqkit stats TSV from HvpAssembly (this.asm_stats_tsv)"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   bin_summary_tsv
        File   asm_stats_tsv
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python3 - <<'PYEOF'
import csv, os, sys
from collections import Counter

bin_file  = "~{bin_summary_tsv}"
asm_file  = "~{asm_stats_tsv}"
out_file  = "~{sample_name}.mag_stats.tsv"
name      = "~{sample_name}"


def gtdb_level(tax_str, prefix):
    if not tax_str or tax_str == "NA":
        return None
    for part in tax_str.split(";"):
        p = part.strip()
        if p.startswith(prefix):
            val = p[len(prefix):]
            return val if val else None
    return None


def mean_r(vals):
    return round(sum(vals) / len(vals), 2) if vals else "NA"


def pct(num, denom):
    try:
        return round(100.0 * num / float(denom), 2) if float(denom) > 0 else "NA"
    except (ValueError, TypeError):
        return "NA"


# Load assembly stats
asm = {}
with open(asm_file) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
    if rows:
        asm = rows[0]

asm_cols = [
    "num_contigs", "bases_in_contigs", "mean_contig_length",
    "q1_contig_length", "median_contig_length", "q3_contig_length",
    "n50_contig_length", "max_contig_length", "mean_contig_gc",
    "num_circ_contigs", "num_1Mb_contigs", "num_circ_1Mb_contigs",
]

# Aggregate bin summary
tier_counts  = Counter()
n_single = n_single_hq = n_no_skani = 0
bases_all = bases_hq = bases_med = 0
comp_all = []; comp_hq = []
cont_all = []; cont_hq = []
qs_hq = []
species_all = set(); species_hq = set()
genera_all  = set(); genera_hq  = set()

with open(bin_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        tier   = row.get("mimag_quality", "NA")
        single = row.get("single_contig", "0") == "1"
        tier_counts[tier] += 1
        if single:
            n_single += 1
            if tier == "high-quality":
                n_single_hq += 1
        if row.get("skani_score", "NA") == "NA":
            n_no_skani += 1
        try:
            gb = int(row.get("genome_size") or 0)
            bases_all += gb
            if tier == "high-quality":
                bases_hq += gb
            elif tier == "medium-quality":
                bases_med += gb
        except ValueError:
            pass
        try:
            comp = float(row["completeness"]); cont = float(row["contamination"])
            qs   = float(row["quality_score"])
            comp_all.append(comp); cont_all.append(cont)
            if tier == "high-quality":
                comp_hq.append(comp); cont_hq.append(cont); qs_hq.append(qs)
        except (ValueError, KeyError):
            pass
        tax = row.get("gtdb_taxonomy", "NA")
        sp  = gtdb_level(tax, "s__"); ge = gtdb_level(tax, "g__")
        if sp:
            species_all.add(sp)
            if tier == "high-quality": species_hq.add(sp)
        if ge:
            genera_all.add(ge)
            if tier == "high-quality": genera_hq.add(ge)

asm_bases  = asm.get("bases_in_contigs")
total_bins = sum(tier_counts.values())

mag_metrics = {
    "total_bins":                total_bins,
    "n_high_quality":            tier_counts["high-quality"],
    "n_medium_quality":          tier_counts["medium-quality"],
    "n_low_quality":             tier_counts["low-quality"],
    "n_contaminated":            tier_counts["contaminated"],
    "n_no_skani_hit":            n_no_skani,
    "n_single_contig":           n_single,
    "n_single_contig_hq":        n_single_hq,
    "total_bases_all_bins":      bases_all,
    "total_bases_hq":            bases_hq,
    "total_bases_medium":        bases_med,
    "pct_assembly_bases_binned": pct(bases_all, asm_bases),
    "pct_assembly_bases_hq":     pct(bases_hq,  asm_bases),
    "mean_completeness_all":     mean_r(comp_all),
    "mean_completeness_hq":      mean_r(comp_hq),
    "mean_contamination_all":    mean_r(cont_all),
    "mean_contamination_hq":     mean_r(cont_hq),
    "mean_quality_score_hq":     mean_r(qs_hq),
    "n_distinct_species_all":    len(species_all),
    "n_distinct_species_hq":     len(species_hq),
    "n_distinct_genera_all":     len(genera_all),
    "n_distinct_genera_hq":      len(genera_hq),
}

all_cols = ["sample_name"] + asm_cols + list(mag_metrics.keys())
with open(out_file, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(all_cols)
    w.writerow([name] + [asm.get(c, "NA") for c in asm_cols] + list(mag_metrics.values()))

# Write scalar txt files for WDL read_int() / read_float()
def stat_f(v): return str(v) if v != "NA" else "0.0"

with open("total_bins.txt",                "w") as f: f.write(str(mag_metrics["total_bins"]))
with open("n_high_quality.txt",            "w") as f: f.write(str(mag_metrics["n_high_quality"]))
with open("n_medium_quality.txt",          "w") as f: f.write(str(mag_metrics["n_medium_quality"]))
with open("n_low_quality.txt",             "w") as f: f.write(str(mag_metrics["n_low_quality"]))
with open("n_contaminated.txt",            "w") as f: f.write(str(mag_metrics["n_contaminated"]))
with open("n_no_skani_hit.txt",            "w") as f: f.write(str(mag_metrics["n_no_skani_hit"]))
with open("n_single_contig.txt",           "w") as f: f.write(str(mag_metrics["n_single_contig"]))
with open("n_single_contig_hq.txt",        "w") as f: f.write(str(mag_metrics["n_single_contig_hq"]))
with open("total_bases_all_bins.txt",      "w") as f: f.write(str(mag_metrics["total_bases_all_bins"]))
with open("total_bases_hq.txt",            "w") as f: f.write(str(mag_metrics["total_bases_hq"]))
with open("total_bases_medium.txt",        "w") as f: f.write(str(mag_metrics["total_bases_medium"]))
with open("pct_assembly_bases_binned.txt", "w") as f: f.write(stat_f(mag_metrics["pct_assembly_bases_binned"]))
with open("pct_assembly_bases_hq.txt",     "w") as f: f.write(stat_f(mag_metrics["pct_assembly_bases_hq"]))
with open("mean_completeness_all.txt",     "w") as f: f.write(stat_f(mag_metrics["mean_completeness_all"]))
with open("mean_completeness_hq.txt",      "w") as f: f.write(stat_f(mag_metrics["mean_completeness_hq"]))
with open("mean_contamination_all.txt",    "w") as f: f.write(stat_f(mag_metrics["mean_contamination_all"]))
with open("mean_contamination_hq.txt",     "w") as f: f.write(stat_f(mag_metrics["mean_contamination_hq"]))
with open("mean_quality_score_hq.txt",     "w") as f: f.write(stat_f(mag_metrics["mean_quality_score_hq"]))
with open("n_distinct_species_all.txt",    "w") as f: f.write(str(mag_metrics["n_distinct_species_all"]))
with open("n_distinct_species_hq.txt",     "w") as f: f.write(str(mag_metrics["n_distinct_species_hq"]))
with open("n_distinct_genera_all.txt",     "w") as f: f.write(str(mag_metrics["n_distinct_genera_all"]))
with open("n_distinct_genera_hq.txt",      "w") as f: f.write(str(mag_metrics["n_distinct_genera_hq"]))

print(f"Wrote {out_file}: {total_bins} bins, "
      f"{tier_counts['high-quality']} HQ, "
      f"{len(species_hq)} distinct HQ species", file=sys.stderr)
PYEOF
    >>>

    output {
        File    mag_stats_tsv             = "~{sample_name}.mag_stats.tsv"
        Int     total_bins                = read_int("total_bins.txt")
        Int     n_high_quality            = read_int("n_high_quality.txt")
        Int     n_medium_quality          = read_int("n_medium_quality.txt")
        Int     n_low_quality             = read_int("n_low_quality.txt")
        Int     n_contaminated            = read_int("n_contaminated.txt")
        Int     n_no_skani_hit            = read_int("n_no_skani_hit.txt")
        Int     n_single_contig           = read_int("n_single_contig.txt")
        Int     n_single_contig_hq        = read_int("n_single_contig_hq.txt")
        Int     total_bases_all_bins      = read_int("total_bases_all_bins.txt")
        Int     total_bases_hq            = read_int("total_bases_hq.txt")
        Int     total_bases_medium        = read_int("total_bases_medium.txt")
        Float   pct_assembly_bases_binned = read_float("pct_assembly_bases_binned.txt")
        Float   pct_assembly_bases_hq     = read_float("pct_assembly_bases_hq.txt")
        Float   mean_completeness_all     = read_float("mean_completeness_all.txt")
        Float   mean_completeness_hq      = read_float("mean_completeness_hq.txt")
        Float   mean_contamination_all    = read_float("mean_contamination_all.txt")
        Float   mean_contamination_hq     = read_float("mean_contamination_hq.txt")
        Float   mean_quality_score_hq     = read_float("mean_quality_score_hq.txt")
        Int     n_distinct_species_all    = read_int("n_distinct_species_all.txt")
        Int     n_distinct_species_hq     = read_int("n_distinct_species_hq.txt")
        Int     n_distinct_genera_all     = read_int("n_distinct_genera_all.txt")
        Int     n_distinct_genera_hq      = read_int("n_distinct_genera_hq.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            4,
        disk_gb:           20,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:          select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:       select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb,  default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:   select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:      select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:       select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:           select_first([runtime_attr.docker,            default_attr.docker])
    }
}
