version 1.0

import "../../structs/Structs.wdl"

task Minimap2AlignReads {

    meta {
        description: "Align PacBio HiFi reads to assembly contigs with minimap2 (map-hifi preset), sort and index the output BAM with samtools, and compute a per-sample alignment summary. Emits the BAM, its index, a flagstat report, a summary TSV, and per-metric scalar outputs."

        tool:          "minimap2 + samtools"
        tool_version:  "minimap2 2.31 / samtools 1.21"
        tool_url:      "https://github.com/lh3/minimap2"
        tool_citation: "Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100."

        outputs: {
            sorted_bam:      "Read-to-contig alignment BAM (minimap2 map-hifi, coordinate-sorted)",
            sorted_bam_bai:  "BAI index for sorted_bam",
            flagstat_txt:    "samtools flagstat output",
            align_stats_tsv: "Alignment summary TSV (total_reads, mapped_reads, fraction_mapped)",
            total_reads:     "Total primary reads in the alignment (from flagstat)",
            mapped_reads:    "Total primary reads mapped to assembly contigs",
            fraction_mapped: "Fraction of primary reads mapped (mapped_reads / total_reads)"
        }
    }

    parameter_meta {
        input_fastq:           "HiFi reads in FASTQ format (gzipped)"
        assembly_fasta:        "Assembly contigs FASTA to align reads against"
        sample_name:           "Sample identifier used as the output file prefix (e.g. bc2097)"
        extra_args:            "Additional command-line args appended verbatim to the minimap2 invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   input_fastq
        File   assembly_fasta
        String sample_name

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_fastq, "GB") + 5.0 * size(assembly_fasta, "GB"))

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

        minimap2 -ax map-hifi -t "${NUM_CPUS}" ~{extra_args} ~{assembly_fasta} ~{input_fastq} \
            | samtools sort -@ 8 -m "${MEM_PER_THREAD_GB}G" -o ~{sample_name}.sorted.bam -

        samtools index ~{sample_name}.sorted.bam

        samtools flagstat ~{sample_name}.sorted.bam > ~{sample_name}.flagstat.txt

        awk '
        /^[0-9]+ \+ [0-9]+ primary$/           { total = $1 }
        /^[0-9]+ \+ [0-9]+ primary mapped /    { mapped = $1 }
        END {
            frac = (total > 0) ? mapped / total : 0
            print "total_reads\tmapped_reads\tfraction_mapped"
            printf "%d\t%d\t%.6f\n", total, mapped, frac
        }' ~{sample_name}.flagstat.txt > ~{sample_name}.align_stats.tsv

        # Extract per-metric scalars from the alignment summary TSV
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["total_reads"])     > "stat.total_reads.txt"
            print $(h["mapped_reads"])    > "stat.mapped_reads.txt"
            print $(h["fraction_mapped"]) > "stat.fraction_mapped.txt"
        }' ~{sample_name}.align_stats.tsv
    >>>

    output {
        File  sorted_bam      = "~{sample_name}.sorted.bam"
        File  sorted_bam_bai  = "~{sample_name}.sorted.bam.bai"
        File  flagstat_txt    = "~{sample_name}.flagstat.txt"
        File  align_stats_tsv = "~{sample_name}.align_stats.tsv"

        Int   total_reads     = read_int("stat.total_reads.txt")
        Int   mapped_reads    = read_int("stat.mapped_reads.txt")
        Float fraction_mapped = read_float("stat.fraction_mapped.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
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
