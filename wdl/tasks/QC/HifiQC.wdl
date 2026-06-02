version 1.0

import "../../structs/Structs.wdl"

task HifiSeqkitStats {

    meta {
        description: "Run seqkit stats and fx2tab on a HiFi FASTQ file to collect read metrics including length distribution, quality thresholds, and bases in reads above 10kb / 20kb cutoffs. Emits both file artifacts and per-metric scalar outputs."

        tool:          "seqkit"
        tool_version:  "2.12.0"
        tool_url:      "https://github.com/shenwei356/seqkit"
        tool_citation: "Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE. 2016;11(10):e0163962."

        outputs: {
            seqkit_stats_tsv:         "Output of seqkit stats -T -a (length stats, N50, Q20/Q30 percentages, GC%)",
            seqkit_fx2tab_tsv:        "Two-row TSV with summed bases in reads >10kb and >20kb",
            num_reads:                "Total number of HiFi reads (seqkit num_seqs)",
            bases_in_reads:           "Total bases across all reads (seqkit sum_len)",
            max_read_length:          "Maximum read length",
            q1_read_length:           "Q1 (25th percentile) of read length distribution",
            median_read_length:       "Median (Q2) read length",
            q3_read_length:           "Q3 (75th percentile) of read length distribution",
            n50_read_length:          "N50 read length",
            pct_q20_bases:            "Percentage of bases at Q20 or higher",
            pct_q30_bases:            "Percentage of bases at Q30 or higher",
            mean_read_gc:             "Mean GC content (percentage) across reads, formatted to one decimal place",
            bases_in_reads_over_10kb: "Sum of read lengths for reads longer than 10kb",
            bases_in_reads_over_20kb: "Sum of read lengths for reads longer than 20kb"
        }
    }

    parameter_meta {
        input_fastq:           "HiFi reads in FASTQ format (gzipped)"
        extra_args:            "Additional command-line args appended verbatim to the seqkit stats invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_fastq

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

        # seqkit stats: key read metrics in tabular format
        seqkit stats -T -a -j "${NUM_CPUS}" ~{extra_args} ~{input_fastq} \
            | cut -f 4,5,7,8,9,10,11,13,15,16,18 \
            > ~{fq_basename}.seqkit_stats.tsv

        # seqkit fx2tab: bases in reads >10kb and >20kb
        seqkit fx2tab -n --length ~{input_fastq} | awk '
        BEGIN {
            OFS = "\t"
            sum_10kb = 0
            sum_20kb = 0
        }
        {
            if ($NF > 10000) {
                sum_10kb += $NF
            }
            if ($NF > 20000) {
                sum_20kb += $NF
            }
        }
        END {
            print "bases_in_reads_over_10kb", "bases_in_reads_over_20kb"
            print sum_10kb, sum_20kb
        }' > ~{fq_basename}.seqkit_fx2tab.tsv

        # Extract per-metric scalars from the seqkit stats TSV (canonical key names + type casts).
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["num_seqs"])         > "stat.num_reads.txt"
            print $(h["sum_len"])          > "stat.bases_in_reads.txt"
            print $(h["max_len"])          > "stat.max_read_length.txt"
            print int($(h["Q1"]))          > "stat.q1_read_length.txt"
            print int($(h["Q2"]))          > "stat.median_read_length.txt"
            print int($(h["Q3"]))          > "stat.q3_read_length.txt"
            print $(h["N50"])              > "stat.n50_read_length.txt"
            print $(h["Q20(%)"])           > "stat.pct_q20_bases.txt"
            print $(h["Q30(%)"])           > "stat.pct_q30_bases.txt"
            printf "%.1f\n", $(h["GC(%)"]) > "stat.mean_read_gc.txt"
        }' ~{fq_basename}.seqkit_stats.tsv

        # Extract per-metric scalars from the fx2tab TSV.
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["bases_in_reads_over_10kb"]) > "stat.bases_in_reads_over_10kb.txt"
            print $(h["bases_in_reads_over_20kb"]) > "stat.bases_in_reads_over_20kb.txt"
        }' ~{fq_basename}.seqkit_fx2tab.tsv
    >>>

    output {
        File seqkit_stats_tsv         = "~{fq_basename}.seqkit_stats.tsv"
        File seqkit_fx2tab_tsv        = "~{fq_basename}.seqkit_fx2tab.tsv"

        Int   num_reads                = read_int("stat.num_reads.txt")
        Int   bases_in_reads           = read_int("stat.bases_in_reads.txt")
        Int   max_read_length          = read_int("stat.max_read_length.txt")
        Int   q1_read_length           = read_int("stat.q1_read_length.txt")
        Int   median_read_length       = read_int("stat.median_read_length.txt")
        Int   q3_read_length           = read_int("stat.q3_read_length.txt")
        Int   n50_read_length          = read_int("stat.n50_read_length.txt")
        Float pct_q20_bases            = read_float("stat.pct_q20_bases.txt")
        Float pct_q30_bases            = read_float("stat.pct_q30_bases.txt")
        Float mean_read_gc             = read_float("stat.mean_read_gc.txt")
        Int   bases_in_reads_over_10kb = read_int("stat.bases_in_reads_over_10kb.txt")
        Int   bases_in_reads_over_20kb = read_int("stat.bases_in_reads_over_20kb.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
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

task HifiKraken2 {

    meta {
        description: "Run kraken2 on HiFi reads and parse the resulting report to extract a classification summary (percent bacteria/virus/fungi/human/unclassified, plus top-ranked genus and species). Emits both file artifacts and per-metric scalar outputs."

        tool:          "kraken2"
        tool_version:  "2.17.1"
        tool_url:      "https://github.com/DerrickWood/kraken2"
        tool_citation: "Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biology. 2019;20(1):257."

        outputs: {
            kraken_report:    "Kraken2 report (--report output) listing per-taxon read counts and percentages",
            kraken_output:    "Per-read kraken2 classification (--output)",
            kraken2_stats:    "Parsed summary TSV with bacteria/virus/fungi/human/unclassified counts and top genus/species",
            pct_bacteria:     "Percentage of reads classified as Bacteria",
            pct_virus:        "Percentage of reads classified as Viruses",
            pct_fungi:        "Percentage of reads classified as Fungi",
            pct_human:        "Percentage of reads classified as Homo sapiens",
            pct_unclassified: "Percentage of reads left unclassified",
            top_genus:        "Genus with the highest kraken2 classification percentage",
            pct_top_genus:    "Percentage of reads assigned to top_genus",
            top_species:      "Species with the highest kraken2 classification percentage",
            pct_top_species:  "Percentage of reads assigned to top_species"
        }
    }

    parameter_meta {
        input_fastq:           "HiFi reads in FASTQ format (gzipped)"
        kraken2_db_hash:       "Kraken2 database hash.k2d file"
        kraken2_db_opts:       "Kraken2 database opts.k2d file"
        kraken2_db_taxo:       "Kraken2 database taxo.k2d file"
        confidence:            "Kraken2 confidence threshold (default 0.001)"
        extra_args:            "Additional command-line args appended verbatim to the kraken2 invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File  input_fastq
        File  kraken2_db_hash
        File  kraken2_db_opts
        File  kraken2_db_taxo
        Float confidence = 0.001

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    String fq_basename = basename(input_fastq, ".fastq.gz")

    # Disk: DB files are read-only and don't transform — keep them at 1x.
    # Only the FASTQ scales by the §6 default 5x multiplier.
    Int disk_size = 20 + ceil(size(kraken2_db_hash, "GB") + size(kraken2_db_opts, "GB") + size(kraken2_db_taxo, "GB")) + ceil(5.0 * size(input_fastq, "GB"))

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

        # Assemble kraken2 database directory from pre-extracted files
        mkdir -p kraken2_db
        ln -s ~{kraken2_db_hash} kraken2_db/hash.k2d
        ln -s ~{kraken2_db_opts} kraken2_db/opts.k2d
        ln -s ~{kraken2_db_taxo} kraken2_db/taxo.k2d

        kraken2 \
            --confidence ~{confidence} \
            --output ~{fq_basename}.kraken_output.txt \
            --report ~{fq_basename}.kraken_report.txt \
            --threads "${NUM_CPUS}" \
            --db kraken2_db \
            ~{extra_args} \
            ~{input_fastq}

        # Parse kraken2 report to extract summary statistics
        awk -F'\t' '
        BEGIN {
            OFS = "\t"
            n_bacteria = 0;    pct_bacteria = 0
            n_virus = 0;       pct_virus = 0
            n_fungi = 0;       pct_fungi = 0
            n_human = 0;       pct_human = 0
            n_unclassified = 0; pct_unclassified = 0
            top_genus = "";    top_genus_pct = 0;    top_genus_val = 0
            top_species = "";  top_species_pct = 0;  top_species_val = 0
        }
        {
            pct = $1 + 0
            level = $4
            tax = $6
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", tax)

            if (tax == "Bacteria")          { pct_bacteria = $1;      n_bacteria = $2 }
            else if (tax == "Fungi")        { pct_fungi = $1;         n_fungi = $2 }
            else if (tax == "Homo sapiens") { pct_human = $1;         n_human = $2 }
            else if (tax == "Viruses")      { pct_virus = $1;         n_virus = $2 }
            else if (tax == "unclassified") { pct_unclassified = $1;  n_unclassified = $2 }

            if (level == "G" && pct > top_genus_val) {
                top_genus = tax
                top_genus_pct = $1
                top_genus_val = pct
            }
            if (level == "S" && pct > top_species_val) {
                top_species = tax
                top_species_pct = $1
                top_species_val = pct
            }
        }
        END {
            print "bacteria", "pct_bacteria", "virus", "pct_virus", "fungi", "pct_fungi", \
                  "human", "pct_human", "unclassified", "pct_unclassified", \
                  "top_genus", "pct_top_genus", "top_species", "pct_top_species"
            print n_bacteria, pct_bacteria, n_virus, pct_virus, n_fungi, pct_fungi, \
                  n_human, pct_human, n_unclassified, pct_unclassified, \
                  top_genus, top_genus_pct, top_species, top_species_pct
        }' ~{fq_basename}.kraken_report.txt > ~{fq_basename}.kraken2_stats.tsv

        # Extract per-metric scalars from the parsed summary TSV.
        awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i }
        NR == 2 {
            print $(h["pct_bacteria"])     > "stat.pct_bacteria.txt"
            print $(h["pct_virus"])        > "stat.pct_virus.txt"
            print $(h["pct_fungi"])        > "stat.pct_fungi.txt"
            print $(h["pct_human"])        > "stat.pct_human.txt"
            print $(h["pct_unclassified"]) > "stat.pct_unclassified.txt"
            print $(h["top_genus"])        > "stat.top_genus.txt"
            print $(h["pct_top_genus"])    > "stat.pct_top_genus.txt"
            print $(h["top_species"])      > "stat.top_species.txt"
            print $(h["pct_top_species"])  > "stat.pct_top_species.txt"
        }' ~{fq_basename}.kraken2_stats.tsv
    >>>

    output {
        File kraken_report = "~{fq_basename}.kraken_report.txt"
        File kraken_output = "~{fq_basename}.kraken_output.txt"
        File kraken2_stats = "~{fq_basename}.kraken2_stats.tsv"

        Float  pct_bacteria     = read_float("stat.pct_bacteria.txt")
        Float  pct_virus        = read_float("stat.pct_virus.txt")
        Float  pct_fungi        = read_float("stat.pct_fungi.txt")
        Float  pct_human        = read_float("stat.pct_human.txt")
        Float  pct_unclassified = read_float("stat.pct_unclassified.txt")
        String top_genus        = read_string("stat.top_genus.txt")
        Float  pct_top_genus    = read_float("stat.pct_top_genus.txt")
        String top_species      = read_string("stat.top_species.txt")
        Float  pct_top_species  = read_float("stat.pct_top_species.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             128,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "staphb/kraken2:2.17.1"
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

task HifiReadStats {

    meta {
        description: "Aggregate pre-computed per-sample scalar metrics into 3 coverage estimates and a human-readable boxed ASCII report. All non-coverage values are passed through from upstream tasks; this task does no parsing of upstream TSVs."

        outputs: {
            hifi_read_stats_report:   "Formatted ASCII-table report grouped by metric category",
            estimate_cvg:             "Estimated coverage = bases_in_reads / expected_genome_size; '-' if expected_genome_size is not available",
            estimate_cvg_reads_10kb:  "Estimated coverage from reads >10kb; '-' if expected_genome_size is not available",
            estimate_cvg_reads_20kb:  "Estimated coverage from reads >20kb; '-' if expected_genome_size is not available"
        }
    }

    parameter_meta {
        sample_name:              "Sample name used as output prefix (e.g. bc2019)"
        genus:                    "Genus name (passthrough into the report)"
        species:                  "Species name (passthrough into the report)"
        tax_id:                   "NCBI taxonomy ID"
        expected_genome_size:     "Expected genome size in bases, or 'NA' if unknown. Numeric values drive coverage estimation."
        num_reads:                "Total number of reads"
        bases_in_reads:           "Total bases across all reads (used for estimate_cvg)"
        bases_in_reads_over_10kb: "Sum of read lengths for reads >10kb (used for estimate_cvg_reads_10kb)"
        bases_in_reads_over_20kb: "Sum of read lengths for reads >20kb (used for estimate_cvg_reads_20kb)"
        q1_read_length:           "Q1 of read length distribution"
        median_read_length:       "Median read length"
        q3_read_length:           "Q3 of read length distribution"
        n50_read_length:          "N50 read length"
        max_read_length:          "Maximum read length"
        pct_q20_bases:            "Percentage of bases at Q20 or higher"
        pct_q30_bases:            "Percentage of bases at Q30 or higher"
        mean_read_accuracy:       "Mean per-read accuracy (percentage)"
        mean_qual_score:          "Mean per-read Phred quality score"
        mean_passes:              "Mean number of CCS subreads per read"
        mean_read_gc:             "Mean GC content (percentage) across reads"
        pct_bacteria:             "Percentage of reads classified as Bacteria"
        pct_fungi:                "Percentage of reads classified as Fungi"
        pct_virus:                "Percentage of reads classified as Viruses"
        pct_human:                "Percentage of reads classified as Homo sapiens"
        pct_unclassified:         "Percentage of reads left unclassified"
        top_genus:                "Genus with the highest kraken2 classification percentage"
        pct_top_genus:            "Percentage of reads assigned to top_genus"
        top_species:              "Species with the highest kraken2 classification percentage"
        pct_top_species:          "Percentage of reads assigned to top_species"
        runtime_attr_override:    "Override the default runtime attributes"
    }

    input {
        String sample_name
        String genus
        String species
        Int    tax_id
        String expected_genome_size

        Int    num_reads
        Int    bases_in_reads
        Int    bases_in_reads_over_10kb
        Int    bases_in_reads_over_20kb
        Int    q1_read_length
        Int    median_read_length
        Int    q3_read_length
        Int    n50_read_length
        Int    max_read_length
        Float  pct_q20_bases
        Float  pct_q30_bases
        Float  mean_read_gc

        Float  mean_read_accuracy
        Float  mean_qual_score
        Int    mean_passes

        Float  pct_bacteria
        Float  pct_fungi
        Float  pct_virus
        Float  pct_human
        Float  pct_unclassified
        String top_genus
        Float  pct_top_genus
        String top_species
        Float  pct_top_species

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

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

        EXPECTED_GENOME_SIZE="~{expected_genome_size}"
        SUM_LEN=~{bases_in_reads}
        B10K=~{bases_in_reads_over_10kb}
        B20K=~{bases_in_reads_over_20kb}

        if [[ "${EXPECTED_GENOME_SIZE}" =~ ^[0-9]+$ ]] && [[ "${EXPECTED_GENOME_SIZE}" -gt 0 ]]; then
            EST_CVG=$(( SUM_LEN / EXPECTED_GENOME_SIZE ))
            EST_CVG_10KB=$(( B10K / EXPECTED_GENOME_SIZE ))
            EST_CVG_20KB=$(( B20K / EXPECTED_GENOME_SIZE ))
        else
            EST_CVG="-"
            EST_CVG_10KB="-"
            EST_CVG_20KB="-"
        fi

        printf '%s\n' "${EST_CVG}"      > estimate_cvg.txt
        printf '%s\n' "${EST_CVG_10KB}" > estimate_cvg_reads_10kb.txt
        printf '%s\n' "${EST_CVG_20KB}" > estimate_cvg_reads_20kb.txt

        REPORT="~{sample_name}.hifi_read_stats.report.txt"
        BAR="+----------------------------+---------------------+"
        FMT="| %-26s | %-19s |\n"

        {
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "Metric" "Value"
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "id" "~{sample_name}"
            # shellcheck disable=SC2059
            printf "${FMT}" "genus" "~{genus}"
            # shellcheck disable=SC2059
            printf "${FMT}" "species" "~{species}"
            # shellcheck disable=SC2059
            printf "${FMT}" "tax_id" "~{tax_id}"
            # shellcheck disable=SC2059
            printf "${FMT}" "expected_genome_size" "${EXPECTED_GENOME_SIZE}"
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "num_reads" "~{num_reads}"
            # shellcheck disable=SC2059
            printf "${FMT}" "bases_in_reads" "${SUM_LEN}"
            # shellcheck disable=SC2059
            printf "${FMT}" "estimate_cvg" "${EST_CVG}"
            # shellcheck disable=SC2059
            printf "${FMT}" "bases_in_reads_over_10kb" "${B10K}"
            # shellcheck disable=SC2059
            printf "${FMT}" "estimate_cvg_reads_10kb" "${EST_CVG_10KB}"
            # shellcheck disable=SC2059
            printf "${FMT}" "bases_in_reads_over_20kb" "${B20K}"
            # shellcheck disable=SC2059
            printf "${FMT}" "estimate_cvg_reads_20kb" "${EST_CVG_20KB}"
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "q1_read_length" "~{q1_read_length}"
            # shellcheck disable=SC2059
            printf "${FMT}" "median_read_length" "~{median_read_length}"
            # shellcheck disable=SC2059
            printf "${FMT}" "q3_read_length" "~{q3_read_length}"
            # shellcheck disable=SC2059
            printf "${FMT}" "n50_read_length" "~{n50_read_length}"
            # shellcheck disable=SC2059
            printf "${FMT}" "max_read_length" "~{max_read_length}"
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_q20_bases" "~{pct_q20_bases}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_q30_bases" "~{pct_q30_bases}"
            # shellcheck disable=SC2059
            printf "${FMT}" "mean_read_accuracy" "~{mean_read_accuracy}"
            # shellcheck disable=SC2059
            printf "${FMT}" "mean_qual_score" "~{mean_qual_score}"
            # shellcheck disable=SC2059
            printf "${FMT}" "mean_passes" "~{mean_passes}"
            # shellcheck disable=SC2059
            printf "${FMT}" "mean_read_gc" "~{mean_read_gc}"
            echo "${BAR}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_bacteria" "~{pct_bacteria}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_fungi" "~{pct_fungi}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_virus" "~{pct_virus}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_human" "~{pct_human}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_unclassified" "~{pct_unclassified}"
            # shellcheck disable=SC2059
            printf "${FMT}" "top_genus" "~{top_genus}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_top_genus" "~{pct_top_genus}"
            # shellcheck disable=SC2059
            printf "${FMT}" "top_species" "~{top_species}"
            # shellcheck disable=SC2059
            printf "${FMT}" "pct_top_species" "~{pct_top_species}"
            echo "${BAR}"
        } > "${REPORT}"
    >>>

    output {
        File   hifi_read_stats_report  = "~{sample_name}.hifi_read_stats.report.txt"
        String estimate_cvg            = read_string("estimate_cvg.txt")
        String estimate_cvg_reads_10kb = read_string("estimate_cvg_reads_10kb.txt")
        String estimate_cvg_reads_20kb = read_string("estimate_cvg_reads_20kb.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:22.04"
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
