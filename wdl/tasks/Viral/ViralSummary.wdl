version 1.0

import "../../structs/Structs.wdl"

task ViralContigSummary {

    meta {
        description: "Build a per-contig viral summary TSV by outer-joining geNomad and VirSorter2 outputs with their per-tool CheckV quality results. Every contig called viral by at least one tool gets one row; tools that did not call a contig produce NA for that tool's columns. VirSorter2 sub-region calls (||1_1 etc.) are collapsed to one row per base contig. Contig IDs are normalized to strip tool-specific suffixes (||full, ||1_1, |provirus_start_stop) before joining. VIBRANT is not included — it has been removed from the pipeline."

        tool:         "python3"
        tool_version: "3.x"

        outputs: {
            viral_contig_summary_tsv: "Per-contig summary TSV joining geNomad and VirSorter2 scores with their CheckV quality results"
        }
    }

    parameter_meta {
        genomad_summary_tsv:   "geNomad virus_summary.tsv from the Genomad task"
        genomad_checkv_tsv:    "CheckV quality_summary.tsv for geNomad viral sequences"
        vs2_score_tsv:         "VirSorter2 final-viral-score.tsv from the VirSorter2 task"
        vs2_checkv_tsv:        "CheckV quality_summary.tsv for VirSorter2 viral sequences"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   genomad_summary_tsv
        File   genomad_checkv_tsv
        File   vs2_score_tsv
        File   vs2_checkv_tsv
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        python3 - <<'PYEOF'
import csv, re, sys
from collections import defaultdict

def normalize_id(s):
    s = str(s)
    s = re.sub(r'\|\|.*$', '', s)
    s = re.sub(r'\|provirus_\d+_\d+$', '', s)
    s = re.sub(r'_fragment_\d+$', '', s)
    s = re.sub(r'[ _]length=\d+.*$', '', s)
    s = re.sub(r'[ _]mult=\S*', '', s)
    return s.strip()

def read_tsv(path):
    try:
        with open(path, newline='') as f:
            return list(csv.DictReader(f, delimiter='\t'))
    except FileNotFoundError:
        print(f"WARNING: file not found, skipping: {path}", file=sys.stderr)
        return []

def load_genomad(path):
    keep = ['virus_score', 'fdr', 'n_hallmarks', 'topology', 'coordinates', 'taxonomy']
    result = {}
    for row in read_tsv(path):
        cid = normalize_id(row['seq_name'])
        if cid not in result:
            result[cid] = {'genomad_' + k: row.get(k, '') for k in keep}
    return result

CHECKV_QUALITY_RANK = {
    'High-quality': 4, 'Medium-quality': 3, 'Low-quality': 2,
    'Not-determined': 1, '': 0,
}

def load_checkv(path, prefix):
    keep_rename = {
        'checkv_quality':      f'{prefix}_quality',
        'miuvig_quality':      f'{prefix}_miuvig_quality',
        'completeness':        f'{prefix}_completeness',
        'completeness_method': f'{prefix}_completeness_method',
        'contamination':       f'{prefix}_contamination',
        'provirus':            f'{prefix}_provirus',
        'proviral_length':     f'{prefix}_proviral_length',
        'contig_length':       f'{prefix}_contig_length',
        'warnings':            f'{prefix}_warnings',
    }
    best = {}
    for row in read_tsv(path):
        orig = row['contig_id']
        cid = normalize_id(orig)
        is_extract = 1 if re.search(r'\|provirus_\d+_\d+$', orig) else 0
        rank = CHECKV_QUALITY_RANK.get(row.get('checkv_quality', ''), 0)
        if cid not in best or (is_extract, rank) > (best[cid][0], best[cid][1]):
            best[cid] = (is_extract, rank, {new: row.get(orig_col, '') for orig_col, new in keep_rename.items()})
    return {cid: data for cid, (_, _, data) in best.items()}

def load_virsorter(path):
    groups = defaultdict(list)
    for row in read_tsv(path):
        cid = normalize_id(row['seqname'])
        m = re.search(r'\|\|(.+)$', row['seqname'])
        region = m.group(1) if m else 'full'
        try:
            score = float(row.get('max_score') or 0)
        except ValueError:
            score = 0.0
        groups[cid].append((region, score, row))
    result = {}
    for cid, entries in groups.items():
        region_scores = ';'.join(f'{r}:{s:.3f}' for r, s, _ in entries)
        entries.sort(key=lambda x: x[1], reverse=True)
        _, _, best = entries[0]
        result[cid] = {
            'vs2_max_score':         best.get('max_score', ''),
            'vs2_max_score_group':   best.get('max_score_group', ''),
            'vs2_hallmark':          best.get('hallmark', ''),
            'vs2_viral_fraction':    best.get('viral', ''),
            'vs2_cellular_fraction': best.get('cellular', ''),
            'vs2_num_regions':       str(len(entries)),
            'vs2_region_scores':     region_scores,
        }
    return result

OUTPUT_COLS = [
    'contig_id', 'contig_length',
    'genomad_virus_score', 'genomad_fdr', 'genomad_n_hallmarks',
    'genomad_topology', 'genomad_coordinates', 'genomad_taxonomy',
    'genomad_checkv_quality', 'genomad_checkv_miuvig_quality',
    'genomad_checkv_completeness', 'genomad_checkv_completeness_method',
    'genomad_checkv_contamination', 'genomad_checkv_provirus',
    'genomad_checkv_proviral_length', 'genomad_checkv_contig_length',
    'genomad_checkv_warnings',
    'vs2_max_score', 'vs2_max_score_group', 'vs2_hallmark',
    'vs2_viral_fraction', 'vs2_cellular_fraction',
    'vs2_num_regions', 'vs2_region_scores',
    'vs2_checkv_quality', 'vs2_checkv_miuvig_quality',
    'vs2_checkv_completeness', 'vs2_checkv_completeness_method',
    'vs2_checkv_contamination', 'vs2_checkv_provirus',
    'vs2_checkv_proviral_length', 'vs2_checkv_contig_length',
    'vs2_checkv_warnings',
    'n_tools', 'tools_calling',
]

genomad  = load_genomad("~{genomad_summary_tsv}")
gcheckv  = load_checkv("~{genomad_checkv_tsv}", 'genomad_checkv')
vs2      = load_virsorter("~{vs2_score_tsv}")
vscheckv = load_checkv("~{vs2_checkv_tsv}", 'vs2_checkv')

all_ids = sorted(set(genomad) | set(vs2))
tool_presence = [('genomad', genomad), ('vs2', vs2)]

if not all_ids:
    print('WARNING: no viral contigs found from any tool — writing empty output', file=sys.stderr)

with open("~{sample_name}.viral_contig_summary.tsv", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=OUTPUT_COLS, delimiter='\t',
                            extrasaction='ignore', restval='NA')
    writer.writeheader()
    for cid in all_ids:
        row = {'contig_id': cid}
        for d in (genomad, gcheckv, vs2, vscheckv):
            row.update(d.get(cid, {}))
        tools = [t for t, d in tool_presence if cid in d]
        row['n_tools']       = str(len(tools))
        row['tools_calling'] = ','.join(tools)
        lengths = []
        for k in ('genomad_checkv_contig_length', 'vs2_checkv_contig_length'):
            v = row.get(k, '')
            if v and v != 'NA':
                try:
                    lengths.append(int(v))
                except ValueError:
                    pass
        row['contig_length'] = str(max(lengths)) if lengths else 'NA'
        writer.writerow(row)

print(f"Wrote ~{sample_name}.viral_contig_summary.tsv: {len(all_ids)} contigs", file=sys.stderr)
PYEOF
    >>>

    output {
        File viral_contig_summary_tsv = "~{sample_name}.viral_contig_summary.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ViralOverallSummary {

    meta {
        description: "Produce a single-row viral summary TSV and scalar outputs for a sample from the per-contig viral_contig_summary.tsv. Aggregates counts by CheckV quality tier (complete, high-quality, medium-quality, low-quality, not-determined), tool agreement (called by both geNomad and VirSorter2 vs. one tool only), topology (DTR, ITR, provirus), and computes N50 and total viral bases. Scalar outputs are written to the Terra data table."

        tool:         "python3"
        tool_version: "3.x"

        outputs: {
            viral_stats_tsv:        "Single-row TSV with all aggregate viral metrics for the sample",
            n_viral_contigs:        "Total viral contigs called by at least one tool",
            n_hq_viral:             "Number of high-quality viral contigs (CheckV High-quality or Complete)",
            n_complete_viral:       "Number of complete viral genomes (CheckV Complete)",
            n_both_tools:           "Number of contigs called viral by both geNomad and VirSorter2",
            n50_viral:              "N50 length of viral contigs (using CheckV contig lengths)",
            total_viral_bases:      "Total bases across all viral contigs",
            n_1tool:                "Number of contigs called viral by exactly one tool",
            n_genomad:              "Number of contigs called viral by geNomad",
            n_vs2:                  "Number of contigs called viral by VirSorter2",
            n_mq:                   "Number of medium-quality viral contigs (CheckV Medium-quality)",
            n_lq:                   "Number of low-quality viral contigs (CheckV Low-quality)",
            n_nd:                   "Number of viral contigs with CheckV quality Not-determined",
            n_no_checkv:            "Number of viral contigs with no CheckV quality result",
            mean_contig_length:     "Mean contig length across all viral contigs",
            n_dtr:                  "Number of viral contigs with DTR (Direct Terminal Repeats) topology",
            n_itr:                  "Number of viral contigs with ITR (Inverted Terminal Repeats) topology",
            n_provirus:             "Number of viral contigs classified as provirus by geNomad",
            n_no_repeats:           "Number of viral contigs with no terminal repeats",
            mean_completeness:      "Mean CheckV completeness across viral contigs with a completeness estimate",
            mean_contamination:     "Mean CheckV contamination across viral contigs with a contamination estimate"
        }
    }

    parameter_meta {
        viral_contig_summary_tsv: "Per-contig viral summary TSV from ViralContigSummary"
        sample_name:              "Sample identifier used as output file prefix"
        runtime_attr_override:    "Override the default runtime attributes"
    }

    input {
        File   viral_contig_summary_tsv
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        python3 - <<'PYEOF'
import csv, sys
from collections import defaultdict

in_file  = "~{viral_contig_summary_tsv}"
out_file = "~{sample_name}.viral_stats.tsv"

QUALITY_RANK = {
    'Complete': 5, 'High-quality': 4, 'Medium-quality': 3,
    'Low-quality': 2, 'Not-determined': 1, 'NA': 0, '': 0,
}

def best_of(*vals):
    best = ''
    for v in vals:
        if QUALITY_RANK.get(v, 0) > QUALITY_RANK.get(best, 0):
            best = v
    return best if best else 'NA'

def max_numeric(*vals):
    result = None
    for v in vals:
        if v and v != 'NA':
            try:
                f = float(v)
                if result is None or f > result:
                    result = f
            except ValueError:
                pass
    return result

def n50(lengths):
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    target = sum(s) / 2
    cumsum = 0
    for length in s:
        cumsum += length
        if cumsum >= target:
            return length
    return 0

def fmt(v, d=1):
    return f'{v:.{d}f}' if v is not None else 'NA'

n_1tool = n_2tool = 0
n_genomad = n_vs2 = 0
qual_counts = defaultdict(int)
lengths = []
completeness_vals = []
contamination_vals = []
topo_counts = defaultdict(int)

with open(in_file) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

for row in rows:
    n = int(row.get('n_tools') or 0)
    if n == 1: n_1tool += 1
    elif n >= 2: n_2tool += 1

    tools = row.get('tools_calling', '')
    if 'genomad' in tools: n_genomad += 1
    if 'vs2'     in tools: n_vs2     += 1

    bq = best_of(
        row.get('genomad_checkv_quality', ''),
        row.get('vs2_checkv_quality', ''),
    )
    qual_counts[bq] += 1

    # Use proviral_length for proviruses so we don't inflate N50/total_bases
    # with the full host contig. Fall back to vs2 provirus flag if genomad
    # didn't call it, then fall back to full contig_length if neither flagged it.
    is_provirus = (
        row.get('genomad_checkv_provirus', '') == 'Yes' or
        row.get('vs2_checkv_provirus',     '') == 'Yes'
    )
    if is_provirus:
        cl = (row.get('genomad_checkv_proviral_length', '') or
              row.get('vs2_checkv_proviral_length',     ''))
    else:
        cl = row.get('contig_length', '')
    if cl and cl != 'NA':
        try:
            lengths.append(int(cl))
        except ValueError:
            pass

    bc = max_numeric(
        row.get('genomad_checkv_completeness', ''),
        row.get('vs2_checkv_completeness', ''),
    )
    if bc is not None:
        completeness_vals.append(bc)

    best_q_rank = 0
    best_contam = None
    for q_col, c_col in [
        ('genomad_checkv_quality', 'genomad_checkv_contamination'),
        ('vs2_checkv_quality',     'vs2_checkv_contamination'),
    ]:
        q = row.get(q_col, '')
        rank = QUALITY_RANK.get(q, 0)
        if rank > best_q_rank:
            c = row.get(c_col, '')
            if c and c != 'NA':
                try:
                    best_contam = float(c)
                    best_q_rank = rank
                except ValueError:
                    pass
    if best_contam is not None:
        contamination_vals.append(best_contam)

    topo = row.get('genomad_topology', '')
    if topo and topo != 'NA':
        topo_counts[topo] += 1

n_contigs = len(rows)
n_hq = qual_counts.get('High-quality', 0) + qual_counts.get('Complete', 0)
total_bases = sum(lengths)

OUTPUT_COLS = [
    'sample_name',
    'n_contigs', 'n_1tool', 'n_both_tools',
    'n_genomad', 'n_vs2',
    'n_complete', 'n_hq', 'n_mq', 'n_lq', 'n_nd', 'n_no_checkv',
    'total_viral_bases', 'n50_viral', 'mean_contig_length',
    'n_dtr', 'n_itr', 'n_provirus', 'n_no_repeats',
    'mean_completeness', 'mean_contamination',
]

row_out = {
    'sample_name':        "~{sample_name}",
    'n_contigs':          n_contigs,
    'n_1tool':            n_1tool,
    'n_both_tools':       n_2tool,
    'n_genomad':          n_genomad,
    'n_vs2':              n_vs2,
    'n_complete':         qual_counts.get('Complete', 0),
    'n_hq':               n_hq,
    'n_mq':               qual_counts.get('Medium-quality', 0),
    'n_lq':               qual_counts.get('Low-quality', 0),
    'n_nd':               qual_counts.get('Not-determined', 0),
    'n_no_checkv':        qual_counts.get('NA', 0),
    'total_viral_bases':  total_bases,
    'n50_viral':          n50(lengths),
    'mean_contig_length': fmt(total_bases / len(lengths) if lengths else None),
    'n_dtr':              topo_counts.get('DTR', 0),
    'n_itr':              topo_counts.get('ITR', 0),
    'n_provirus':         topo_counts.get('Provirus', 0),
    'n_no_repeats':       topo_counts.get('No terminal repeats', 0),
    'mean_completeness':  fmt(sum(completeness_vals) / len(completeness_vals) if completeness_vals else None),
    'mean_contamination': fmt(sum(contamination_vals) / len(contamination_vals) if contamination_vals else None, 2),
}

with open(out_file, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=OUTPUT_COLS, delimiter='\t')
    w.writeheader()
    w.writerow(row_out)

# Scalar files for WDL read_int / read_float
with open('stat.n_viral_contigs.txt',  'w') as f: f.write(str(n_contigs))
with open('stat.n_hq_viral.txt',       'w') as f: f.write(str(n_hq))
with open('stat.n_complete_viral.txt', 'w') as f: f.write(str(qual_counts.get('Complete', 0)))
with open('stat.n_both_tools.txt',     'w') as f: f.write(str(n_2tool))
with open('stat.n50_viral.txt',        'w') as f: f.write(str(n50(lengths)))
with open('stat.total_viral_bases.txt','w') as f: f.write(str(total_bases))
with open('stat.n_1tool.txt',          'w') as f: f.write(str(n_1tool))
with open('stat.n_genomad.txt',        'w') as f: f.write(str(n_genomad))
with open('stat.n_vs2.txt',            'w') as f: f.write(str(n_vs2))
with open('stat.n_mq.txt',             'w') as f: f.write(str(qual_counts.get('Medium-quality', 0)))
with open('stat.n_lq.txt',             'w') as f: f.write(str(qual_counts.get('Low-quality', 0)))
with open('stat.n_nd.txt',             'w') as f: f.write(str(qual_counts.get('Not-determined', 0)))
with open('stat.n_no_checkv.txt',      'w') as f: f.write(str(qual_counts.get('NA', 0)))
with open('stat.n_dtr.txt',            'w') as f: f.write(str(topo_counts.get('DTR', 0)))
with open('stat.n_itr.txt',            'w') as f: f.write(str(topo_counts.get('ITR', 0)))
with open('stat.n_provirus.txt',       'w') as f: f.write(str(topo_counts.get('Provirus', 0)))
with open('stat.n_no_repeats.txt',     'w') as f: f.write(str(topo_counts.get('No terminal repeats', 0)))
mean_cl = total_bases / len(lengths) if lengths else None
with open('stat.mean_contig_length.txt', 'w') as f: f.write(f'{mean_cl:.1f}' if mean_cl is not None else '0.0')
mean_comp = sum(completeness_vals) / len(completeness_vals) if completeness_vals else None
with open('stat.mean_completeness.txt', 'w') as f: f.write(f'{mean_comp:.1f}' if mean_comp is not None else '0.0')
mean_contam = sum(contamination_vals) / len(contamination_vals) if contamination_vals else None
with open('stat.mean_contamination.txt','w') as f: f.write(f'{mean_contam:.2f}' if mean_contam is not None else '0.0')

print(f"Sample ~{sample_name}: {n_contigs} viral contigs, {n_hq} HQ/complete, "
      f"{n_2tool} called by both tools", file=sys.stderr)
PYEOF
    >>>

    output {
        File  viral_stats_tsv   = "~{sample_name}.viral_stats.tsv"
        Int   n_viral_contigs   = read_int("stat.n_viral_contigs.txt")
        Int   n_hq_viral        = read_int("stat.n_hq_viral.txt")
        Int   n_complete_viral  = read_int("stat.n_complete_viral.txt")
        Int   n_both_tools      = read_int("stat.n_both_tools.txt")
        Int   n50_viral         = read_int("stat.n50_viral.txt")
        Int   total_viral_bases = read_int("stat.total_viral_bases.txt")
        Int   n_1tool           = read_int("stat.n_1tool.txt")
        Int   n_genomad         = read_int("stat.n_genomad.txt")
        Int   n_vs2             = read_int("stat.n_vs2.txt")
        Int   n_mq              = read_int("stat.n_mq.txt")
        Int   n_lq              = read_int("stat.n_lq.txt")
        Int   n_nd              = read_int("stat.n_nd.txt")
        Int   n_no_checkv       = read_int("stat.n_no_checkv.txt")
        Float mean_contig_length = read_float("stat.mean_contig_length.txt")
        Int   n_dtr             = read_int("stat.n_dtr.txt")
        Int   n_itr             = read_int("stat.n_itr.txt")
        Int   n_provirus        = read_int("stat.n_provirus.txt")
        Int   n_no_repeats      = read_int("stat.n_no_repeats.txt")
        Float mean_completeness  = read_float("stat.mean_completeness.txt")
        Float mean_contamination = read_float("stat.mean_contamination.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
