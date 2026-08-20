version 1.0

import "../../structs/Structs.wdl"

task SpikeInRemoval {

    meta {
        description: "Align PacBio HiFi reads to a spike-in reference with minimap2 (map-hifi, --secondary=no), filter hits requiring ≥97% identity and ≥80% query coverage, and remove matching reads from the input FASTQ. Emits a spike-in-cleaned FASTQ, the spike-in read name list (consumed by FilterCleanBam to remove those reads from the original BAM), per-reference count reports, and aggregate scalar metrics. Spike-in references are the NEB HiFi spike-in panel (11 sequences: Lambda phage size ladder, pBR322, phiX174, M13mp18)."

        tool:          "minimap2 + samtools"
        tool_version:  "minimap2 2.31 / samtools 1.21"
        tool_url:      "https://github.com/lh3/minimap2"
        tool_citation: "Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100."

        outputs: {
            cleaned_fastq:                    "Spike-in-removed FASTQ (gzipped); fed into HifiKraken2",
            spikein_read_names:               "One read name per line for reads identified as spike-in; consumed by FilterCleanBam",
            spikein_report:                   "Long-format report: one row per spike-in reference with read count and fraction of total",
            spikein_stats_tsv:                "Wide-format TSV: one column pair (num_<ref>, fraction_<ref>) per spike-in reference",
            total_spikein_reads:              "Total number of unique read names identified as spike-in",
            fraction_spikein:                 "Fraction of input reads identified as spike-in",
            num_Lambda_NEB_2026_125_bp:       "Reads classified as Lambda phage 125bp spike-in",
            fraction_Lambda_NEB_2026_125_bp:  "Fraction of reads classified as Lambda phage 125bp spike-in",
            num_Lambda_NEB_2026_564_bp:       "Reads classified as Lambda phage 564bp spike-in",
            fraction_Lambda_NEB_2026_564_bp:  "Fraction of reads classified as Lambda phage 564bp spike-in",
            num_Lambda_NEB_2026_2027_bp:      "Reads classified as Lambda phage 2027bp spike-in",
            fraction_Lambda_NEB_2026_2027_bp: "Fraction of reads classified as Lambda phage 2027bp spike-in",
            num_Lambda_NEB_2026_2322_bp:      "Reads classified as Lambda phage 2322bp spike-in",
            fraction_Lambda_NEB_2026_2322_bp: "Fraction of reads classified as Lambda phage 2322bp spike-in",
            num_Lambda_NEB_2026_4361_bp:      "Reads classified as Lambda phage 4361bp spike-in",
            fraction_Lambda_NEB_2026_4361_bp: "Fraction of reads classified as Lambda phage 4361bp spike-in",
            num_Lambda_NEB_2026_6557_bp:      "Reads classified as Lambda phage 6557bp spike-in",
            fraction_Lambda_NEB_2026_6557_bp: "Fraction of reads classified as Lambda phage 6557bp spike-in",
            num_Lambda_NEB_2026_9416_bp:      "Reads classified as Lambda phage 9416bp spike-in",
            fraction_Lambda_NEB_2026_9416_bp: "Fraction of reads classified as Lambda phage 9416bp spike-in",
            num_Lambda_NEB_2026_23130_bp:     "Reads classified as Lambda phage 23130bp spike-in",
            fraction_Lambda_NEB_2026_23130_bp:"Fraction of reads classified as Lambda phage 23130bp spike-in",
            num_pBR322_BamHI:                 "Reads classified as pBR322 spike-in",
            fraction_pBR322_BamHI:            "Fraction of reads classified as pBR322 spike-in",
            num_phiX174_NEB_PstI:             "Reads classified as phiX174 spike-in",
            fraction_phiX174_NEB_PstI:        "Fraction of reads classified as phiX174 spike-in",
            num_M13mp18_PstI:                 "Reads classified as M13mp18 spike-in",
            fraction_M13mp18_PstI:            "Fraction of reads classified as M13mp18 spike-in"
        }
    }

    parameter_meta {
        input_fastq:           "HiFi reads in FASTQ format (gzipped)"
        spikein_fasta:         "Spike-in reference FASTA (NEB HiFi spike-in panel, 11 sequences)"
        extra_args:            "Additional command-line args appended verbatim to the minimap2 invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_fastq
        File spikein_fasta

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    String fq_basename = basename(input_fastq, ".fastq.gz")

    Int disk_size = 10 + ceil(5.0 * size(input_fastq, "GB"))

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

        # Align to spike-in reference; --secondary=no retains only the single best hit per read
        minimap2 -ax map-hifi --secondary=no -t "${NUM_CPUS}" ~{extra_args} \
            ~{spikein_fasta} ~{input_fastq} \
            | samtools view -F 4 -h > ~{fq_basename}.spikein.sam

        total_reads=$(zcat ~{input_fastq} | awk 'END{print NR/4}')

        # Parse SAM: seed ref_counts from @SQ headers so all references appear in output.
        # Filter aligned reads requiring ≥97% identity and ≥80% query coverage.
        # Writes matching read names to spikein_readnames.tmp; produces long and wide reports.
        awk -v total="${total_reads}" \
            -v outtxt="~{fq_basename}.spikein_stats.txt" \
            -v outtsv="~{fq_basename}.spikein_stats.tsv" \
            -v outnm="spikein_readnames.tmp" '
        BEGIN { FS="\t" }
        /^@SQ/ {
            for (i=2; i<=NF; i++) if ($i ~ /^SN:/) ref_counts[substr($i,4)] = 0
            next
        }
        /^@/ { next }
        {
            if (and($2, 4)) next
            qname=$1; ref=$3; cigar=$6
            if (cigar == "*") next

            m=0; ins=0; d=0; s=0; h=0
            c=cigar
            while (match(c, /[0-9]+[MIDNSHP=X]/)) {
                seg=substr(c, RSTART, RLENGTH)
                n=substr(seg, 1, length(seg)-1)+0
                op=substr(seg, length(seg), 1)
                if      (op=="M") m  +=n
                else if (op=="I") ins+=n
                else if (op=="D") d  +=n
                else if (op=="S") s  +=n
                else if (op=="H") h  +=n
                c=substr(c, RSTART+RLENGTH)
            }

            nm=0
            for (i=12; i<=NF; i++) if ($i ~ /^NM:i:/) { nm=substr($i,6)+0; break }

            align_len = m+ins+d
            query_len = m+ins+s+h
            if (align_len==0 || query_len==0) next

            identity = (align_len - nm) / align_len
            coverage = (m+ins) / query_len

            if (identity >= 0.97 && coverage >= 0.80) {
                ref_counts[ref]++
                print qname >> outnm
            }
        }
        END {
            n = total+0
            nrefs = split("Lambda_NEB_2026_125_bp Lambda_NEB_2026_564_bp Lambda_NEB_2026_2027_bp Lambda_NEB_2026_2322_bp Lambda_NEB_2026_4361_bp Lambda_NEB_2026_6557_bp Lambda_NEB_2026_9416_bp Lambda_NEB_2026_23130_bp pBR322_BamHI phiX174_NEB_PstI M13mp18_PstI", refs, " ")

            print "spikein_reference\tread_count\tfraction_of_total" > outtxt
            header = ""; vals = ""
            for (i=1; i<=nrefs; i++) {
                r = refs[i]
                cnt = (r in ref_counts) ? ref_counts[r] : 0
                frac = (n > 0) ? cnt / n : 0
                printf "%s\t%d\t%.6f\n", r, cnt, frac >> outtxt
                printf "%d\n",   cnt  > ("stat.num_" r ".txt")
                printf "%.6f\n", frac > ("stat.fraction_" r ".txt")
                sep = (i > 1) ? "\t" : ""
                header = header sep "num_" r "\tfraction_" r
                vals   = vals   sep cnt "\t" sprintf("%.6f", frac)
            }
            print header > outtsv
            print vals >> outtsv
        }' ~{fq_basename}.spikein.sam

        rm ~{fq_basename}.spikein.sam

        # Deduplicate read names (a read may appear in multiple supplementary alignments)
        if [[ -f spikein_readnames.tmp ]]; then
            sort -u spikein_readnames.tmp > ~{fq_basename}.spikein_readnames.txt
            rm spikein_readnames.tmp
        else
            touch ~{fq_basename}.spikein_readnames.txt
        fi

        # Aggregate scalar metrics
        awk 'END{print NR}' ~{fq_basename}.spikein_readnames.txt > stat.total_spikein_reads.txt
        awk -v n="${total_reads}" \
            'END{ printf "%.6f\n", (n > 0) ? NR / n : 0 }' \
            ~{fq_basename}.spikein_readnames.txt > stat.fraction_spikein.txt

        # Remove spike-in reads from FASTQ; skip filtering pass if none were found
        if [[ ! -s ~{fq_basename}.spikein_readnames.txt ]]; then
            cp ~{input_fastq} ~{fq_basename}.spikein_cleaned.fastq.gz
        else
            zcat ~{input_fastq} | awk -v namefile="~{fq_basename}.spikein_readnames.txt" '
            BEGIN {
                while ((getline line < namefile) > 0) names[line] = 1
                close(namefile)
            }
            FNR%4==1 {
                split(substr($0,2), parts, " ")
                skip = (parts[1] in names)
            }
            !skip { print }
            ' | gzip > ~{fq_basename}.spikein_cleaned.fastq.gz
        fi
    >>>

    output {
        File  cleaned_fastq       = "~{fq_basename}.spikein_cleaned.fastq.gz"
        File  spikein_read_names  = "~{fq_basename}.spikein_readnames.txt"
        File  spikein_report      = "~{fq_basename}.spikein_stats.txt"
        File  spikein_stats_tsv   = "~{fq_basename}.spikein_stats.tsv"
        Int   total_spikein_reads = read_int("stat.total_spikein_reads.txt")
        Float fraction_spikein    = read_float("stat.fraction_spikein.txt")

        Int   num_Lambda_NEB_2026_125_bp       = read_int("stat.num_Lambda_NEB_2026_125_bp.txt")
        Float fraction_Lambda_NEB_2026_125_bp  = read_float("stat.fraction_Lambda_NEB_2026_125_bp.txt")
        Int   num_Lambda_NEB_2026_564_bp       = read_int("stat.num_Lambda_NEB_2026_564_bp.txt")
        Float fraction_Lambda_NEB_2026_564_bp  = read_float("stat.fraction_Lambda_NEB_2026_564_bp.txt")
        Int   num_Lambda_NEB_2026_2027_bp      = read_int("stat.num_Lambda_NEB_2026_2027_bp.txt")
        Float fraction_Lambda_NEB_2026_2027_bp = read_float("stat.fraction_Lambda_NEB_2026_2027_bp.txt")
        Int   num_Lambda_NEB_2026_2322_bp      = read_int("stat.num_Lambda_NEB_2026_2322_bp.txt")
        Float fraction_Lambda_NEB_2026_2322_bp = read_float("stat.fraction_Lambda_NEB_2026_2322_bp.txt")
        Int   num_Lambda_NEB_2026_4361_bp      = read_int("stat.num_Lambda_NEB_2026_4361_bp.txt")
        Float fraction_Lambda_NEB_2026_4361_bp = read_float("stat.fraction_Lambda_NEB_2026_4361_bp.txt")
        Int   num_Lambda_NEB_2026_6557_bp      = read_int("stat.num_Lambda_NEB_2026_6557_bp.txt")
        Float fraction_Lambda_NEB_2026_6557_bp = read_float("stat.fraction_Lambda_NEB_2026_6557_bp.txt")
        Int   num_Lambda_NEB_2026_9416_bp      = read_int("stat.num_Lambda_NEB_2026_9416_bp.txt")
        Float fraction_Lambda_NEB_2026_9416_bp = read_float("stat.fraction_Lambda_NEB_2026_9416_bp.txt")
        Int   num_Lambda_NEB_2026_23130_bp     = read_int("stat.num_Lambda_NEB_2026_23130_bp.txt")
        Float fraction_Lambda_NEB_2026_23130_bp = read_float("stat.fraction_Lambda_NEB_2026_23130_bp.txt")
        Int   num_pBR322_BamHI                 = read_int("stat.num_pBR322_BamHI.txt")
        Float fraction_pBR322_BamHI            = read_float("stat.fraction_pBR322_BamHI.txt")
        Int   num_phiX174_NEB_PstI             = read_int("stat.num_phiX174_NEB_PstI.txt")
        Float fraction_phiX174_NEB_PstI        = read_float("stat.fraction_phiX174_NEB_PstI.txt")
        Int   num_M13mp18_PstI                 = read_int("stat.num_M13mp18_PstI.txt")
        Float fraction_M13mp18_PstI            = read_float("stat.fraction_M13mp18_PstI.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-align:0.1.0"
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

task FilterCleanBam {

    meta {
        description: "Remove spike-in and human reads from the original unmapped BAM by read name exclusion. Spike-in read names come from SpikeInRemoval; human read names (taxid 9606) are extracted here from the kraken2 per-read output. Filtering is applied directly to the original BAM so all PacBio-specific tags (base modifications, rq:f, np:i) are preserved in the output. Produces an indexed cleaned BAM that is the final output of HvpReadProcessing and the primary input to HvpAssembly."

        tool:          "samtools"
        tool_version:  "1.21"
        tool_url:      "https://www.htslib.org/"
        tool_citation: "Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008."

        outputs: {
            cleaned_bam:           "BAM with spike-in and human reads removed; all original BAM tags preserved",
            cleaned_bam_bai:       "BAI index for cleaned_bam",
            num_spikein_excluded:  "Number of unique read names excluded as spike-in",
            num_human_excluded:    "Number of unique read names excluded as human (taxid 9606)",
            num_reads_cleaned_bam: "Total reads remaining in cleaned_bam (samtools view -c)"
        }
    }

    parameter_meta {
        input_bam:             "Original unmapped PacBio HiFi BAM (before any read filtering)"
        spikein_read_names:    "Spike-in read name list produced by SpikeInRemoval (one name per line)"
        kraken_output:         "Per-read kraken2 classification file from HifiKraken2"
        sample_name:           "Sample identifier used as output file prefix (e.g. bc2097)"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   input_bam
        File   spikein_read_names
        File   kraken_output
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bam, "GB"))

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

        # Extract human read names (taxid 9606) from the kraken2 per-read output
        awk -F'\t' '$1=="C" && $3==9606 { print $2 }' ~{kraken_output} \
            | sort -u > human_readnames.txt

        # Count excluded reads per category
        awk 'END{print NR}' ~{spikein_read_names} > stat.num_spikein_excluded.txt
        awk 'END{print NR}' human_readnames.txt    > stat.num_human_excluded.txt

        # Merge into a single exclusion set
        cat ~{spikein_read_names} human_readnames.txt | sort -u > exclude_readnames.txt

        # Stream BAM through awk to exclude reads by name, then re-encode as BAM
        samtools view -h ~{input_bam} \
            | awk 'BEGIN { while ((getline < "exclude_readnames.txt") > 0) excl[$0] = 1 }
                   /^@/ { print; next }
                   !($1 in excl) { print }' \
            | samtools view -b -@ "${NUM_CPUS}" -o ~{sample_name}.cleaned.bam -

        samtools index ~{sample_name}.cleaned.bam
        samtools view -c ~{sample_name}.cleaned.bam > stat.num_reads_cleaned_bam.txt
    >>>

    output {
        File cleaned_bam           = "~{sample_name}.cleaned.bam"
        File cleaned_bam_bai       = "~{sample_name}.cleaned.bam.bai"
        Int  num_spikein_excluded  = read_int("stat.num_spikein_excluded.txt")
        Int  num_human_excluded    = read_int("stat.num_human_excluded.txt")
        Int  num_reads_cleaned_bam = read_int("stat.num_reads_cleaned_bam.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "quay.io/biocontainers/samtools:1.23--h96c455f_0"
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
