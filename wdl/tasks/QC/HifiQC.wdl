version 1.0

import "../../structs/Structs.wdl"

task HifiSeqkitStats {

    meta {
        description: "Run seqkit stats and fx2tab on a HiFi FASTQ file to collect read metrics including length distribution, quality thresholds, and bases in reads above 10kb / 20kb cutoffs."

        tool:          "seqkit"
        tool_version:  "2.12.0"
        tool_url:      "https://github.com/shenwei356/seqkit"
        tool_citation: "Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE. 2016;11(10):e0163962."

        outputs: {
            seqkit_stats_tsv:  "Output of seqkit stats -T -a (length stats, N50, Q20/Q30 percentages, GC%)",
            seqkit_fx2tab_tsv: "Two-row TSV with summed bases in reads >10kb and >20kb"
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
    >>>

    output {
        File seqkit_stats_tsv  = "~{fq_basename}.seqkit_stats.tsv"
        File seqkit_fx2tab_tsv = "~{fq_basename}.seqkit_fx2tab.tsv"
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
        description: "Run kraken2 on HiFi reads and parse the resulting report to extract a classification summary (percent bacteria/virus/fungi/human/unclassified, plus top-ranked genus and species)."

        tool:          "kraken2"
        tool_version:  "2.17.1"
        tool_url:      "https://github.com/DerrickWood/kraken2"
        tool_citation: "Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biology. 2019;20(1):257."

        outputs: {
            kraken_report:  "Kraken2 report (--report output) listing per-taxon read counts and percentages",
            kraken_output:  "Per-read kraken2 classification (--output)",
            kraken2_stats:  "Parsed summary TSV with bacteria/virus/fungi/human/unclassified counts and top genus/species"
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
    >>>

    output {
        File kraken_report = "~{fq_basename}.kraken_report.txt"
        File kraken_output = "~{fq_basename}.kraken_output.txt"
        File kraken2_stats = "~{fq_basename}.kraken2_stats.tsv"
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
        description: "Compile per-sample summary metrics from the bam, seqkit, kraken2, and taxonomy outputs into a human-readable formatted report plus individual scalar outputs for each metric."

        outputs: {
            hifi_read_stats_report:   "Formatted ASCII-table report grouped by metric category",
            num_reads:                "Total number of HiFi reads (from seqkit num_seqs)",
            bases_in_reads:           "Total bases across all reads (from seqkit sum_len)",
            estimate_cvg:             "Estimated coverage = bases_in_reads / expected_genome_size; '-' if expected_genome_size is not available",
            bases_in_reads_over_10kb: "Sum of read lengths for reads longer than 10kb",
            estimate_cvg_reads_10kb:  "Estimated coverage from reads >10kb; '-' if expected_genome_size is not available",
            bases_in_reads_over_20kb: "Sum of read lengths for reads longer than 20kb",
            estimate_cvg_reads_20kb:  "Estimated coverage from reads >20kb; '-' if expected_genome_size is not available",
            q1_read_length:           "Q1 (25th percentile) of read length distribution",
            median_read_length:       "Median (Q2) read length",
            q3_read_length:           "Q3 (75th percentile) of read length distribution",
            n50_read_length:          "N50 read length",
            max_read_length:          "Maximum read length",
            pct_q20_bases:            "Percentage of bases at Q20 or higher",
            pct_q30_bases:            "Percentage of bases at Q30 or higher",
            mean_read_base_qual:      "Mean per-read Phred quality score derived from rq:f tags",
            mean_read_accuracy:       "Mean per-read accuracy (as a percentage) derived from rq:f tags",
            mean_num_passes:          "Mean number of CCS subreads per read (from np:i tags)",
            mean_read_gc:             "Mean GC content (percentage) across reads",
            pct_bacteria:             "Percentage of reads classified as Bacteria by kraken2",
            pct_fungi:                "Percentage of reads classified as Fungi by kraken2",
            pct_virus:                "Percentage of reads classified as Viruses by kraken2",
            pct_human:                "Percentage of reads classified as Homo sapiens by kraken2",
            pct_unclassified:         "Percentage of reads left unclassified by kraken2",
            top_genus:                "Genus with the highest kraken2 classification percentage",
            pct_top_genus:            "Percentage of reads assigned to top_genus",
            top_species:              "Species with the highest kraken2 classification percentage",
            pct_top_species:          "Percentage of reads assigned to top_species"
        }
    }

    parameter_meta {
        sample_name:               "Sample name used as output prefix (e.g. bc2019)"
        bam_stats_tsv:             "Output from BamToFastqAndStats task"
        seqkit_stats_tsv:          "Output from HifiSeqkitStats task (stats)"
        seqkit_fx2tab_tsv:         "Output from HifiSeqkitStats task (fx2tab)"
        kraken2_stats_tsv:         "Output from HifiKraken2 task"
        taxid_and_genome_size_tsv: "Output from GetTaxIdAndGenomeSize task"
        extra_args:                "Additional command-line args appended verbatim to the paste invocation"
        runtime_attr_override:     "Override the default runtime attributes"
    }

    input {
        String sample_name
        File   bam_stats_tsv
        File   seqkit_stats_tsv
        File   seqkit_fx2tab_tsv
        File   kraken2_stats_tsv
        File   taxid_and_genome_size_tsv

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(bam_stats_tsv, "GB") + size(seqkit_stats_tsv, "GB") + size(seqkit_fx2tab_tsv, "GB") + size(kraken2_stats_tsv, "GB") + size(taxid_and_genome_size_tsv, "GB")))

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

        paste ~{extra_args} ~{seqkit_stats_tsv} ~{seqkit_fx2tab_tsv} ~{bam_stats_tsv} ~{kraken2_stats_tsv} ~{taxid_and_genome_size_tsv} | awk -v id="~{sample_name}" '
        BEGIN {
            FS = "\t"
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "num_seqs") num_seqs_col = i
                if ($i == "sum_len") sum_len_col = i
                if ($i == "avg_len") avg_len_col = i
                if ($i == "max_len") max_len_col = i
                if ($i == "Q1") q1_col = i
                if ($i == "Q2") q2_col = i
                if ($i == "Q3") q3_col = i
                if ($i == "N50") n50_col = i
                if ($i == "Q20(%)") q20_col = i
                if ($i == "Q30(%)") q30_col = i
                if ($i == "GC(%)") gc_col = i
                if ($i == "bases_in_reads_over_10kb") bases_10kb_col = i
                if ($i == "bases_in_reads_over_20kb") bases_20kb_col = i
                if ($i == "mean_read_accuracy") read_accuracy_col = i
                if ($i == "mean_qual_score") qual_score_col = i
                if ($i == "mean_passes") passes_col = i
                if ($i == "pct_bacteria") pct_bacteria_col = i
                if ($i == "pct_virus") pct_virus_col = i
                if ($i == "pct_fungi") pct_fungi_col = i
                if ($i == "pct_human") pct_human_col = i
                if ($i == "pct_unclassified") pct_unclassified_col = i
                if ($i == "top_genus") top_genus_col = i
                if ($i == "pct_top_genus") pct_top_genus_col = i
                if ($i == "top_species") top_species_col = i
                if ($i == "pct_top_species") pct_top_species_col = i
                if ($i == "tax_id") tax_id_col = i
                if ($i == "expected_genome_size") expected_genome_size_col = i
                if ($i == "genus") genus_col = i
                if ($i == "species") species_col = i
            }
        }
        NR == 2 {
            num_seqs = $(num_seqs_col)
            sum_len = $(sum_len_col)
            avg_len = $(avg_len_col)
            max_len = $(max_len_col)
            q1 = int($(q1_col))
            median_len = int($(q2_col))
            q3 = int($(q3_col))
            n50 = $(n50_col)
            q20_pct = $(q20_col)
            q30_pct = $(q30_col)
            avg_gc = sprintf("%.1f", $(gc_col))
            bases_in_reads_over_10kb = $(bases_10kb_col)
            bases_in_reads_over_20kb = $(bases_20kb_col)
            mean_read_accuracy = $(read_accuracy_col)
            mean_qual_score = $(qual_score_col)
            mean_passes = $(passes_col)
            tax_id = $(tax_id_col)
            expected_genome_size = $(expected_genome_size_col)
            genus = $(genus_col)
            species = $(species_col)

            if (expected_genome_size ~ /^[0-9]/ && expected_genome_size + 0 > 0) {
                est_cvg = int(sum_len / expected_genome_size)
                est_cvg_10kb = (bases_in_reads_over_10kb + 0 > 0) ? int(bases_in_reads_over_10kb / expected_genome_size) : 0
                est_cvg_20kb = (bases_in_reads_over_20kb + 0 > 0) ? int(bases_in_reads_over_20kb / expected_genome_size) : 0
            } else {
                est_cvg = "-"
                est_cvg_10kb = "-"
                est_cvg_20kb = "-"
            }

            pct_bacteria = $(pct_bacteria_col)
            pct_virus = $(pct_virus_col)
            pct_fungi = $(pct_fungi_col)
            pct_human = $(pct_human_col)
            pct_unclassified = $(pct_unclassified_col)
            top_genus = $(top_genus_col)
            pct_top_genus = $(pct_top_genus_col)
            top_species = $(top_species_col)
            pct_top_species = $(pct_top_species_col)
        }
        END {
            out_txt = id ".hifi_read_stats.txt"
            print "id:", id > out_txt
            print "genus:", genus >> out_txt
            print "species:", species >> out_txt
            print "tax_id:", tax_id >> out_txt
            print "expected_genome_size:", expected_genome_size >> out_txt
            print "num_reads:", num_seqs >> out_txt
            print "bases_in_reads:", sum_len >> out_txt
            print "estimate_cvg:", est_cvg >> out_txt
            print "bases_in_reads_over_10kb:", bases_in_reads_over_10kb >> out_txt
            print "estimate_cvg_reads_10kb:", est_cvg_10kb >> out_txt
            print "bases_in_reads_over_20kb:", bases_in_reads_over_20kb >> out_txt
            print "estimate_cvg_reads_20kb:", est_cvg_20kb >> out_txt
            print "q1_read_length:", q1 >> out_txt
            print "median_read_length:", median_len >> out_txt
            print "q3_read_length:", q3 >> out_txt
            print "n50_read_length:", n50 >> out_txt
            print "max_read_length:", max_len >> out_txt
            print "pct_q20_bases:", q20_pct >> out_txt
            print "pct_q30_bases:", q30_pct >> out_txt
            print "mean_read_base_qual:", mean_qual_score >> out_txt
            print "mean_read_accuracy:", mean_read_accuracy >> out_txt
            print "mean_num_passes:", mean_passes >> out_txt
            print "mean_read_gc:", avg_gc >> out_txt
            print "pct_bacteria:", pct_bacteria >> out_txt
            print "pct_fungi:", pct_fungi >> out_txt
            print "pct_virus:", pct_virus >> out_txt
            print "pct_human:", pct_human >> out_txt
            print "pct_unclassified:", pct_unclassified >> out_txt
            print "top_genus:", top_genus >> out_txt
            print "pct_top_genus:", pct_top_genus >> out_txt
            print "top_species:", top_species >> out_txt
            print "pct_top_species:", pct_top_species >> out_txt

            out_tsv = id ".hifi_read_stats.tsv"
            print "id", "genus", "species", "tax_id", "expected_genome_size", \
                  "num_reads", "bases_in_reads", "estimate_cvg", \
                  "bases_in_reads_over_10kb", "estimate_cvg_reads_10kb", \
                  "bases_in_reads_over_20kb", "estimate_cvg_reads_20kb", \
                  "q1_read_length", "median_read_length", "q3_read_length", \
                  "n50_read_length", "max_read_length", \
                  "pct_q20_bases", "pct_q30_bases", "mean_read_base_qual", \
                  "mean_read_accuracy", "mean_num_passes", "mean_read_gc", \
                  "pct_bacteria", "pct_fungi", "pct_virus", "pct_human", "pct_unclassified", \
                  "top_genus", "pct_top_genus", "top_species", "pct_top_species" > out_tsv

            print id, genus, species, tax_id, expected_genome_size, \
                  num_seqs, sum_len, est_cvg, \
                  bases_in_reads_over_10kb, est_cvg_10kb, \
                  bases_in_reads_over_20kb, est_cvg_20kb, \
                  q1, median_len, q3, n50, max_len, \
                  q20_pct, q30_pct, mean_qual_score, \
                  mean_read_accuracy, mean_passes, avg_gc, \
                  pct_bacteria, pct_fungi, pct_virus, pct_human, pct_unclassified, \
                  top_genus, pct_top_genus, top_species, pct_top_species >> out_tsv

            out_report = id ".hifi_read_stats.report.txt"
            bar = "+----------------------------+---------------------+"
            fmt = "| %-26s | %-19s |\n"
            print bar > out_report
            printf fmt, "Metric", "Value" >> out_report
            print bar >> out_report
            printf fmt, "id", id >> out_report
            printf fmt, "genus", genus >> out_report
            printf fmt, "species", species >> out_report
            printf fmt, "tax_id", tax_id >> out_report
            printf fmt, "expected_genome_size", expected_genome_size >> out_report
            print bar >> out_report
            printf fmt, "num_reads", num_seqs >> out_report
            printf fmt, "bases_in_reads", sum_len >> out_report
            printf fmt, "estimate_cvg", est_cvg >> out_report
            printf fmt, "bases_in_reads_over_10kb", bases_in_reads_over_10kb >> out_report
            printf fmt, "estimate_cvg_reads_10kb", est_cvg_10kb >> out_report
            printf fmt, "bases_in_reads_over_20kb", bases_in_reads_over_20kb >> out_report
            printf fmt, "estimate_cvg_reads_20kb", est_cvg_20kb >> out_report
            print bar >> out_report
            printf fmt, "q1_read_length", q1 >> out_report
            printf fmt, "median_read_length", median_len >> out_report
            printf fmt, "q3_read_length", q3 >> out_report
            printf fmt, "n50_read_length", n50 >> out_report
            printf fmt, "max_read_length", max_len >> out_report
            print bar >> out_report
            printf fmt, "pct_q20_bases", q20_pct >> out_report
            printf fmt, "pct_q30_bases", q30_pct >> out_report
            printf fmt, "mean_read_base_qual", mean_qual_score >> out_report
            printf fmt, "mean_read_accuracy", mean_read_accuracy >> out_report
            printf fmt, "mean_num_passes", mean_passes >> out_report
            printf fmt, "mean_read_gc", avg_gc >> out_report
            print bar >> out_report
            printf fmt, "pct_bacteria", pct_bacteria >> out_report
            printf fmt, "pct_fungi", pct_fungi >> out_report
            printf fmt, "pct_virus", pct_virus >> out_report
            printf fmt, "pct_human", pct_human >> out_report
            printf fmt, "pct_unclassified", pct_unclassified >> out_report
            printf fmt, "top_genus", top_genus >> out_report
            printf fmt, "pct_top_genus", pct_top_genus >> out_report
            printf fmt, "top_species", top_species >> out_report
            printf fmt, "pct_top_species", pct_top_species >> out_report
            print bar >> out_report
        }'

        # Split the single-row TSV into per-column scalar files for WDL output bindings.
        tsv="~{sample_name}.hifi_read_stats.tsv"
        header=$(head -1 "$tsv")
        data=$(tail -1 "$tsv")
        IFS=$'\t' read -ra H <<< "$header"
        IFS=$'\t' read -ra D <<< "$data"
        for i in "${!H[@]}"; do
            printf '%s\n' "${D[$i]}" > "stat.${H[$i]}.txt"
        done
    >>>

    output {
        File   hifi_read_stats_report   = "~{sample_name}.hifi_read_stats.report.txt"

        Int    num_reads                = read_int("stat.num_reads.txt")
        Int    bases_in_reads           = read_int("stat.bases_in_reads.txt")
        String estimate_cvg             = read_string("stat.estimate_cvg.txt")
        Int    bases_in_reads_over_10kb = read_int("stat.bases_in_reads_over_10kb.txt")
        String estimate_cvg_reads_10kb  = read_string("stat.estimate_cvg_reads_10kb.txt")
        Int    bases_in_reads_over_20kb = read_int("stat.bases_in_reads_over_20kb.txt")
        String estimate_cvg_reads_20kb  = read_string("stat.estimate_cvg_reads_20kb.txt")
        Int    q1_read_length           = read_int("stat.q1_read_length.txt")
        Int    median_read_length       = read_int("stat.median_read_length.txt")
        Int    q3_read_length           = read_int("stat.q3_read_length.txt")
        Int    n50_read_length          = read_int("stat.n50_read_length.txt")
        Int    max_read_length          = read_int("stat.max_read_length.txt")
        Float  pct_q20_bases            = read_float("stat.pct_q20_bases.txt")
        Float  pct_q30_bases            = read_float("stat.pct_q30_bases.txt")
        Float  mean_read_base_qual      = read_float("stat.mean_read_base_qual.txt")
        Float  mean_read_accuracy       = read_float("stat.mean_read_accuracy.txt")
        Int    mean_num_passes          = read_int("stat.mean_num_passes.txt")
        Float  mean_read_gc             = read_float("stat.mean_read_gc.txt")
        Float  pct_bacteria             = read_float("stat.pct_bacteria.txt")
        Float  pct_fungi                = read_float("stat.pct_fungi.txt")
        Float  pct_virus                = read_float("stat.pct_virus.txt")
        Float  pct_human                = read_float("stat.pct_human.txt")
        Float  pct_unclassified         = read_float("stat.pct_unclassified.txt")
        String top_genus                = read_string("stat.top_genus.txt")
        Float  pct_top_genus            = read_float("stat.pct_top_genus.txt")
        String top_species              = read_string("stat.top_species.txt")
        Float  pct_top_species          = read_float("stat.pct_top_species.txt")
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
