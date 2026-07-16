version 1.0

import "../../structs/Structs.wdl"

task ReadRescue {

    meta {
        description: "Identify reads 'unaccounted for' by the assembly — unmapped or with (soft+hard)-clip fraction exceeding clip_frac_max — and partition them by kraken2 classification into viral and unclassified rescue sets. Cross-references per-read kraken2 output against the minimap2 read-to-contig BAM. The cleaned BAM is converted to FASTQ internally, and matching reads are extracted as gzipped FASTA for downstream viral discovery tools (geNomad, VirSorter2). MAPQ is not used as an unaccounted criterion: low MAPQ indicates multi-mapping, not poor alignment quality."

        tool:         "samtools + seqkit + python3"
        tool_version: "samtools 1.21 / seqkit 2.12.0 / python3"
        tool_url:     "https://www.htslib.org/"

        outputs: {
            inventory_tsv:              "Per-read TSV: read_id, kraken_status, kraken_taxid, kraken_taxon_name, mapped, mapq, clip_frac, read_len, accounted, category, contig",
            rescue_viral_fa_gz:         "Gzipped FASTA of rescued kraken-viral reads (unaccounted and length >= min_rescue_read_len)",
            rescue_unclassified_fa_gz:  "Gzipped FASTA of rescued kraken-unclassified reads (unaccounted and length >= min_rescue_read_len)",
            rescue_summary_tsv:         "Summary TSV: per-category read counts and threshold parameters used",
            num_rescue_viral:           "Number of reads rescued as kraken2-viral",
            num_rescue_unclassified:    "Number of reads rescued as kraken2-unclassified"
        }
    }

    parameter_meta {
        sorted_bam:            "Read-to-contig alignment BAM from HvpAssembly (minimap2 map-hifi, coordinate-sorted); used for clip fraction analysis"
        sorted_bam_bai:        "BAI index for sorted_bam"
        cleaned_bam:           "Cleaned unmapped BAM from HvpReadProcessing; converted to FASTQ internally for rescue read extraction"
        kraken_output:         "Per-read kraken2 classification file from HvpReadProcessing (tab-separated: status, read_id, taxid, ...)"
        kraken_report:         "Kraken2 report from HvpReadProcessing; defines the viral taxid subtree (reads the Viruses clade and all descendants)"
        sample_name:           "Sample identifier used as output file prefix (e.g. bc2097)"
        clip_frac_max:         "Maximum (soft+hard) clip fraction for a read to be considered accounted-for by the assembly (default 0.5)"
        min_rescue_read_len:   "Minimum read length in bases for a rescued read to be included in the output FASTA (default 1000)"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   sorted_bam
        File   sorted_bam_bai
        File   cleaned_bam
        File   kraken_output
        File   kraken_report
        String sample_name

        Float  clip_frac_max       = 0.5
        Int    min_rescue_read_len = 1000

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(cleaned_bam, "GB")) + ceil(2.0 * size(sorted_bam, "GB"))

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        # Convert cleaned BAM → FASTQ for seqkit rescue read extraction
        samtools fastq -@ "${NUM_CPUS}" ~{cleaned_bam} | gzip > cleaned.fastq.gz

        # Cross-reference kraken2 classifications with assembly BAM clip fractions
        python3 - <<'PYEOF'
import csv, os, re, subprocess, sys
from collections import Counter

k_out         = "~{kraken_output}"
k_rep         = "~{kraken_report}"
bam           = "~{sorted_bam}"
clip_frac_max = float("~{clip_frac_max}")
min_len       = int("~{min_rescue_read_len}")

# --- 1. Parse kraken report: viral taxid subtree + taxid→name map ---
viral_taxids = set()
taxid_name = {}
with open(k_rep) as f:
    in_viral = False
    viral_indent = None
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        name_raw = parts[5]
        name = name_raw.strip()
        indent = len(name_raw) - len(name_raw.lstrip(" "))
        taxid = parts[4]
        taxid_name[taxid] = name
        if not in_viral:
            if name == "Viruses":
                in_viral = True
                viral_indent = indent
                viral_taxids.add(taxid)
        else:
            if indent > viral_indent:
                viral_taxids.add(taxid)
            else:
                in_viral = False
print(f"viral taxids in report: {len(viral_taxids)}", file=sys.stderr)

# --- 2. Per-read kraken status ---
read_kraken = {}
with open(k_out) as f:
    for line in f:
        c, rid, taxid, *_ = line.rstrip("\n").split("\t")
        if c == "U" or taxid == "0":
            status = "unclassified"
        elif taxid in viral_taxids:
            status = "viral"
        else:
            status = "other"
        read_kraken[rid] = (status, taxid)
print(f"kraken reads: {len(read_kraken)}", file=sys.stderr)

# --- 3. Walk BAM (primary only), compute clip fractions, write inventory + id lists ---
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
def clip_fraction(cigar):
    if cigar == "*":
        return 0.0
    qlen = 0
    clip = 0
    for n, op in CIGAR_OP.findall(cigar):
        n = int(n)
        if op in ("M", "I", "S", "=", "X", "H"):
            qlen += n
        if op in ("S", "H"):
            clip += n
    return clip / qlen if qlen > 0 else 0.0

inv_f = open("inventory.tsv", "w", newline="")
inv_w = csv.writer(inv_f, delimiter="\t")
inv_w.writerow(["read_id", "kraken_status", "kraken_taxid", "kraken_taxon_name",
                "mapped", "mapq", "clip_frac", "read_len", "accounted", "category", "contig"])
vf = open("rescue.viral.ids", "w")
uf = open("rescue.unclassified.ids", "w")
counts = Counter()

p = subprocess.Popen(
    ["samtools", "view", "-F", "0x900", bam],
    stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
)
seen = 0
for line in p.stdout:
    cols = line.split("\t", 11)
    rid      = cols[0]
    flag     = int(cols[1])
    rname    = cols[2]
    mapq     = int(cols[4])
    cigar    = cols[5]
    seq      = cols[9]
    read_len = len(seq) if seq != "*" else 0
    unmapped = bool(flag & 0x4)

    if unmapped:
        clipf, category, accounted = 0.0, "unmapped", False
    else:
        clipf = clip_fraction(cigar)
        if clipf > clip_frac_max:
            category, accounted = "high_clip", False
        else:
            category, accounted = "mapped", True

    kstat, ktax = read_kraken.get(rid, ("other", "0"))
    kname       = taxid_name.get(ktax, "")
    contig_name = "" if unmapped else rname
    inv_w.writerow([rid, kstat, ktax, kname,
                    "N" if unmapped else "Y",
                    mapq, f"{clipf:.3f}", read_len,
                    "Y" if accounted else "N", category, contig_name])
    counts[(kstat, "accounted" if accounted else "unaccounted")] += 1
    counts[(kstat, category)] += 1

    if (not accounted) and read_len >= min_len:
        if kstat == "viral":
            vf.write(rid + "\n")
        elif kstat == "unclassified":
            uf.write(rid + "\n")
    seen += 1

p.wait()
if p.returncode != 0:
    sys.exit(f"samtools view exited {p.returncode}")
inv_f.close(); vf.close(); uf.close()
print(f"bam reads (primary): {seen}", file=sys.stderr)

# --- 4. Summary ---
with open("rescue.summary.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["clip_frac_max", clip_frac_max])
    w.writerow(["min_rescue_read_len", min_len])
    w.writerow(["viral_taxids_in_report", len(viral_taxids)])
    w.writerow(["reads_in_kraken", len(read_kraken)])
    w.writerow(["reads_in_bam_primary", seen])
    for status in ("viral", "unclassified", "other"):
        for tag in ("accounted", "unaccounted", "mapped", "unmapped", "high_clip"):
            w.writerow([f"{status}_{tag}", counts.get((status, tag), 0)])
PYEOF

        # Extract rescue reads from cleaned FASTQ → gzipped FASTA
        if [[ -s rescue.viral.ids ]]; then
            seqkit grep -f rescue.viral.ids cleaned.fastq.gz \
                | awk 'NR%4==1{sub(/^@/,">"); print} NR%4==2{print}' \
                | gzip > ~{sample_name}.rescue_viral.fasta.gz
        else
            : | gzip > ~{sample_name}.rescue_viral.fasta.gz
        fi

        if [[ -s rescue.unclassified.ids ]]; then
            seqkit grep -f rescue.unclassified.ids cleaned.fastq.gz \
                | awk 'NR%4==1{sub(/^@/,">"); print} NR%4==2{print}' \
                | gzip > ~{sample_name}.rescue_unclassified.fasta.gz
        else
            : | gzip > ~{sample_name}.rescue_unclassified.fasta.gz
        fi

        mv inventory.tsv     ~{sample_name}.rescue_inventory.tsv
        mv rescue.summary.tsv ~{sample_name}.rescue_summary.tsv

        awk 'END{print NR}' rescue.viral.ids        > stat.num_rescue_viral.txt
        awk 'END{print NR}' rescue.unclassified.ids > stat.num_rescue_unclassified.txt
    >>>

    output {
        File inventory_tsv             = "~{sample_name}.rescue_inventory.tsv"
        File rescue_viral_fa_gz        = "~{sample_name}.rescue_viral.fasta.gz"
        File rescue_unclassified_fa_gz = "~{sample_name}.rescue_unclassified.fasta.gz"
        File rescue_summary_tsv        = "~{sample_name}.rescue_summary.tsv"
        Int  num_rescue_viral          = read_int("stat.num_rescue_viral.txt")
        Int  num_rescue_unclassified   = read_int("stat.num_rescue_unclassified.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
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
