version 1.0

import "../../structs/Structs.wdl"

task Myloasm {

    meta {
        description: "Assemble metagenomic PacBio HiFi reads with myloasm. Emits a single primary assembly FASTA (assembly_primary.fa)."

        tool:         "myloasm"
        tool_version: "0.5.1"
        tool_url:     "https://github.com/lbcb-sci/myloasm"

        outputs: {
            assembly_primary_fa: "Primary assembly FASTA produced by myloasm (assembly_primary.fa)"
        }
    }

    parameter_meta {
        input_fastq:           "HiFi reads in FASTQ format (gzipped)"
        extra_args:            "Additional command-line args appended verbatim to the myloasm invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_fastq

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(input_fastq, "GB"))

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

        mkdir -p myloasm_out
        myloasm ~{input_fastq} -o myloasm_out -t "${NUM_CPUS}" --hifi ~{extra_args} --clean-dir
    >>>

    output {
        File assembly_primary_fa = "myloasm_out/assembly_primary.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             128,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/myloasm:0.5.1"
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

task SeqkitAssemblyStats {

    meta {
        description: "Filter an assembly FASTA by minimum contig length, then run seqkit stats and fx2tab to compute contig-level quality metrics. Emits both a summary TSV and per-metric scalar outputs."

        tool:          "seqkit"
        tool_version:  "2.12.0"
        tool_url:      "https://github.com/shenwei356/seqkit"
        tool_citation: "Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE. 2016;11(10):e0163962."

        outputs: {
            asm_stats_tsv:         "Wide-format TSV (header row + values row) with all 12 assembly metrics",
            num_contigs:           "Number of contigs passing the min_contig_len filter",
            bases_in_contigs:      "Total bases across all passing contigs",
            mean_contig_length:    "Mean contig length",
            q1_contig_length:      "Q1 (25th percentile) contig length",
            median_contig_length:  "Median contig length",
            q3_contig_length:      "Q3 (75th percentile) contig length",
            n50_contig_length:     "N50 contig length",
            max_contig_length:     "Maximum contig length",
            mean_contig_gc:        "Mean GC content (percentage) across contigs",
            num_circ_contigs:      "Number of circular contigs (myloasm: _circular-yes; metaMDBG: circular=yes; hifiasm_meta: ctg ID ending in 'c')",
            num_1Mb_contigs:       "Number of contigs longer than 1 Mb",
            num_circ_1Mb_contigs:  "Number of contigs that are both circular and longer than 1 Mb"
        }
    }

    parameter_meta {
        assembly_fasta:        "Assembly contigs FASTA (plain or gzipped)"
        min_contig_len:        "Minimum contig length to retain before computing stats; 0 = no filtering"
        extra_args:            "Additional command-line args appended verbatim to the seqkit stats invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File assembly_fasta
        Int  min_contig_len = 0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    String fa_basename = basename(basename(assembly_fasta, ".gz"), ".fa")

    Int disk_size = 10 + ceil(5.0 * size(assembly_fasta, "GB"))

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

        seqkit seq --min-len ~{min_contig_len} ~{assembly_fasta} > filtered_contigs.fasta

        seqkit stats -T -a -j "${NUM_CPUS}" ~{extra_args} filtered_contigs.fasta > stats_raw.tsv

        seqkit fx2tab -n --length filtered_contigs.fasta > fx2tab.tsv

        # Parse seqkit stats columns by header name and rename to output metric names
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        NR == 2 {
            print "num_contigs\tbases_in_contigs\tmean_contig_length\tq1_contig_length\tmedian_contig_length\tq3_contig_length\tn50_contig_length\tmax_contig_length\tmean_contig_gc"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                $col["num_seqs"], $col["sum_len"], $col["avg_len"],
                $col["Q1"], $col["Q2"], $col["Q3"],
                $col["N50"], $col["max_len"], $col["GC(%)"]
        }' stats_raw.tsv > stats_parsed.tsv

        # Extract per-metric scalars from the renamed stats TSV
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["num_contigs"])                    > "stat.num_contigs.txt"
            print $(h["bases_in_contigs"])               > "stat.bases_in_contigs.txt"
            printf "%.2f\n", $(h["mean_contig_length"])  > "stat.mean_contig_length.txt"
            print int($(h["q1_contig_length"]))          > "stat.q1_contig_length.txt"
            print int($(h["median_contig_length"]))      > "stat.median_contig_length.txt"
            print int($(h["q3_contig_length"]))          > "stat.q3_contig_length.txt"
            print $(h["n50_contig_length"])              > "stat.n50_contig_length.txt"
            print $(h["max_contig_length"])              > "stat.max_contig_length.txt"
            printf "%.1f\n", $(h["mean_contig_gc"])      > "stat.mean_contig_gc.txt"
        }' stats_parsed.tsv

        # Count circular and >1 Mb contigs from the fx2tab full-header column
        awk -F'\t' '
        BEGIN { n_circ = 0; n_1mb = 0; n_circ_1mb = 0 }
        {
            circ = ($1 ~ /circular[=\-]yes/ || $1 ~ /\.ctg[0-9]+c$/)
            big  = ($2 + 0 > 1000000)
            if (circ)        n_circ++
            if (big)         n_1mb++
            if (circ && big) n_circ_1mb++
        }
        END {
            print "num_circ_contigs\tnum_1Mb_contigs\tnum_circ_1Mb_contigs"
            printf "%d\t%d\t%d\n", n_circ, n_1mb, n_circ_1mb
        }' fx2tab.tsv > counts.tsv

        # Extract per-metric scalars from the counts TSV
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["num_circ_contigs"])     > "stat.num_circ_contigs.txt"
            print $(h["num_1Mb_contigs"])      > "stat.num_1Mb_contigs.txt"
            print $(h["num_circ_1Mb_contigs"]) > "stat.num_circ_1Mb_contigs.txt"
        }' counts.tsv

        paste stats_parsed.tsv counts.tsv > ~{fa_basename}.asm_stats.tsv
    >>>

    output {
        File  asm_stats_tsv        = "~{fa_basename}.asm_stats.tsv"

        Int   num_contigs          = read_int("stat.num_contigs.txt")
        Int   bases_in_contigs     = read_int("stat.bases_in_contigs.txt")
        Float mean_contig_length   = read_float("stat.mean_contig_length.txt")
        Int   q1_contig_length     = read_int("stat.q1_contig_length.txt")
        Int   median_contig_length = read_int("stat.median_contig_length.txt")
        Int   q3_contig_length     = read_int("stat.q3_contig_length.txt")
        Int   n50_contig_length    = read_int("stat.n50_contig_length.txt")
        Int   max_contig_length    = read_int("stat.max_contig_length.txt")
        Float mean_contig_gc       = read_float("stat.mean_contig_gc.txt")
        Int   num_circ_contigs     = read_int("stat.num_circ_contigs.txt")
        Int   num_1Mb_contigs      = read_int("stat.num_1Mb_contigs.txt")
        Int   num_circ_1Mb_contigs = read_int("stat.num_circ_1Mb_contigs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "staphb/seqkit:2.12.0"
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
