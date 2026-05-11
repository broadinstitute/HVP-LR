version 1.0

# ============================================================================
# Tasks
# ============================================================================

task parse_sample_sheet {
    meta {
        description: "Parse a sample sheet TSV (with header) and extract columns as parallel arrays."
    }

    parameter_meta {
        # note: for metagenomic samples we will not have genus, species. Instead just the BAM file and some sample identifier
        sample_tsv: {
            description: "Input TSV with columns: bam, genus, species, strain"
        }
    }

    input {
        File sample_tsv
    }

    command <<<
        set -euo pipefail
        tail -n +2 ~{sample_tsv} | cut -f1 > bams.txt
        tail -n +2 ~{sample_tsv} | cut -f2 > genera.txt
        tail -n +2 ~{sample_tsv} | cut -f3 > species.txt
        tail -n +2 ~{sample_tsv} | cut -f4 > strains.txt
        tail -n +2 ~{sample_tsv} | cut -f1 | xargs -I{} basename {} .bam | grep -o 'bc[0-9]*' > barcodes.txt
    >>>

    output {
        Array[String] bam_paths = read_lines("bams.txt")
        Array[String] genera    = read_lines("genera.txt")
        Array[String] species   = read_lines("species.txt")
        Array[String] strains   = read_lines("strains.txt")
        Array[String] barcodes  = read_lines("barcodes.txt")
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        disks:  "local-disk 1 SSD"
        cpu:    1
    }
}

task bam_to_fastq_and_stats {
    meta {
        description: "Convert a BAM file to fastq.gz and calculate mean read accuracy, Phred quality score, and number of passes."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
    }

    parameter_meta {
        input_bam: {
            description: "PacBio HiFi reads BAM file",
            patterns: ["*.bam"]
        }
    }

    input {
        File input_bam

        String  docker_image = "quay.io/biocontainers/samtools:1.23--h96c455f_0"
        Int     disk_size_gb = ceil(3 * size(input_bam, "GB")) + 10
        Int     memory_gb    = 4
        Int     cpu          = 4
    }

    String bam_basename = basename(input_bam, ".bam")

    command <<<
        set -euo pipefail

        samtools fastq \
            -@ ~{cpu - 1} \
            ~{input_bam} \
            | gzip > ~{bam_basename}.fastq.gz

        samtools view ~{input_bam} | awk '
        BEGIN {
            FS = "\t"
            OFS = "\t"
            total_passes = 0
            total_accuracy = 0
            read_count = 0
        }
        {
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^np:i:/) {
                    total_passes += substr($i, 6)
                } else if ($i ~ /^rq:f:/) {
                    total_accuracy += substr($i, 6)
                }
            }
            read_count++
        }
        END {
            if (read_count > 0) {
                mean_passes = int(total_passes / read_count)
                mean_read_accuracy = total_accuracy / read_count
                p_error = 1 - mean_read_accuracy
                mean_phred_score = sprintf("%.1f", -10 * (log(p_error) / log(10)))
                mean_read_accuracy = sprintf("%.3f", 100 * mean_read_accuracy)
            } else {
                mean_read_accuracy = 0
                mean_phred_score = 0
                mean_passes = 0
            }
            print "mean_read_accuracy", "mean_qual_score", "mean_passes"
            print mean_read_accuracy, mean_phred_score, mean_passes
        }' > ~{bam_basename}.bam_stats.tsv

        # Write individual outputs for downstream use
        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f1 > mean_read_accuracy.txt
        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f2 > mean_qual_score.txt
        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f3 > mean_passes.txt
    >>>

    output {
        File   fastq_gz           = "~{bam_basename}.fastq.gz"
        File   bam_stats_tsv      = "~{bam_basename}.bam_stats.tsv"
        Float  mean_read_accuracy = read_float("mean_read_accuracy.txt")
        Float  mean_qual_score    = read_float("mean_qual_score.txt")
        Int    mean_passes        = read_int("mean_passes.txt")
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks:  "local-disk ~{disk_size_gb} SSD"
        cpu:    cpu
    }
}

task hifi_seqkit_stats {
    meta {
        description: "Run seqkit stats and fx2tab on a HiFi FASTQ file to collect read metrics including length distribution and base thresholds."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
    }

    parameter_meta {
        input_fastq: {
            description: "HiFi reads in FASTQ format (gzipped)",
            patterns: ["*.fastq.gz", "*.fq.gz"]
        }
    }

    input {
        File input_fastq

        String  docker_image = "staphb/seqkit:2.12.0"
        Int     disk_size_gb = ceil(2 * size(input_fastq, "GB")) + 10
        Int     memory_gb    = 4
        Int     cpu          = 2
    }

    String fq_basename = basename(input_fastq, ".fastq.gz")

    command <<<
        set -euo pipefail

        # seqkit stats: key read metrics in tabular format
        seqkit stats -T -a ~{input_fastq} \
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

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks:  "local-disk ~{disk_size_gb} SSD"
        cpu:    cpu
    }
}

task hifi_kraken2 {
    meta {
        description: "Run kraken2 on HiFi reads and parse the report to extract classification summary statistics."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
    }

    parameter_meta {
        input_fastq: {
            description: "HiFi reads in FASTQ format (gzipped)",
            patterns: ["*.fastq.gz", "*.fq.gz"]
        }
        kraken2_db_hash: {
            description: "Kraken2 database hash.k2d file"
        }
        kraken2_db_opts: {
            description: "Kraken2 database opts.k2d file"
        }
        kraken2_db_taxo: {
            description: "Kraken2 database taxo.k2d file"
        }
    }

    input {
        File   input_fastq
        File   kraken2_db_hash
        File   kraken2_db_opts
        File   kraken2_db_taxo
        Float  confidence  = 0.001

        String  docker_image = "staphb/kraken2:2.17.1"
        Int     disk_size_gb = ceil(size(input_fastq, "GB") + size(kraken2_db_hash, "GB") + size(kraken2_db_opts, "GB") + size(kraken2_db_taxo, "GB")) + 20
        Int     memory_gb    = 128
        Int     cpu          = 8
    }

    String fq_basename = basename(input_fastq, ".fastq.gz")

    command <<<
        set -euo pipefail

        # Assemble kraken2 database directory from pre-extracted files
        mkdir -p kraken2_db
        ln -s ~{kraken2_db_hash} kraken2_db/hash.k2d
        ln -s ~{kraken2_db_opts} kraken2_db/opts.k2d
        ln -s ~{kraken2_db_taxo} kraken2_db/taxo.k2d

        # Run kraken2
        kraken2 \
            --confidence ~{confidence} \
            --output ~{fq_basename}.kraken_output.txt \
            --report ~{fq_basename}.kraken_report.txt \
            --threads ~{cpu} \
            --db kraken2_db \
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
        File kraken_report    = "~{fq_basename}.kraken_report.txt"
        File kraken_output    = "~{fq_basename}.kraken_output.txt"
        File kraken2_stats    = "~{fq_basename}.kraken2_stats.tsv"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks:  "local-disk ~{disk_size_gb} SSD"
        cpu:    cpu
    }
}

# note: skip for metagenomic samples
task batch_get_taxID_and_genome_size {
    meta {
        description: "Look up NCBI taxonomy IDs and expected genome sizes for all samples in a single VM, deduplicating by genus+species."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
    }

    parameter_meta {
        genera: {
            description: "Array of genus names, one per sample"
        }
        species_list: {
            description: "Array of species names, one per sample"
        }
        barcodes: {
            description: "Array of sample barcodes, one per sample"
        }
    }

    input {
        Array[String] genera
        Array[String] species_list
        Array[String] barcodes

        String  docker_image = "quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1"
        Int     disk_size_gb = 10
        Int     memory_gb    = 2
        Int     cpu          = 1
    }

    File genera_file   = write_lines(genera)
    File species_file  = write_lines(species_list)
    File barcodes_file = write_lines(barcodes)

    command <<<
        set -euo pipefail

        # Download and extract NCBI taxonomy database (once)
        mkdir -p taxdump
        wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf taxdump.tar.gz -C taxdump

        # Download species genome size table (once)
        wget -qO species_genome_size.txt \
            "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt"

        mkdir -p output_tsvs

        # Associative arrays to cache lookups by genus+species
        declare -A taxid_cache
        declare -A genome_size_cache

        # Process each sample
        paste ~{genera_file} ~{species_file} ~{barcodes_file} | while IFS=$'\t' read -r genus species barcode; do
            species_name="${genus} ${species}"
            cache_key="${genus}_${species}"

            # Look up tax_id (cached per unique genus+species)
            if [[ -z "${taxid_cache[$cache_key]+x}" ]]; then
                tax_id=$(echo "$species_name" \
                    | taxonkit name2taxid --sci-name --data-dir taxdump \
                    | cut -f2)

                if [ -z "$tax_id" ]; then
                    tax_id=$(echo "$species_name" \
                        | taxonkit name2taxid --data-dir taxdump \
                        | cut -f2)
                fi

                if [ -z "$tax_id" ]; then
                    echo "Error: no taxonomy ID found for '$species_name'" >&2
                    exit 1
                fi

                taxid_cache[$cache_key]="$tax_id"

                # Look up genome size
                expected_genome_size=$(awk -F'\t' -v id="$tax_id" '$1 == id {print $4}' species_genome_size.txt)
                if [ -z "$expected_genome_size" ]; then
                    echo "Warning: no genome size entry found for tax_id $tax_id. Setting to NA." >&2
                    expected_genome_size="NA"
                fi
                genome_size_cache[$cache_key]="$expected_genome_size"
            fi

            tax_id="${taxid_cache[$cache_key]}"
            expected_genome_size="${genome_size_cache[$cache_key]}"

            # Write per-sample TSV
            outfile="output_tsvs/${barcode}.taxid_and_genome_size.tsv"
            printf "tax_id\texpected_genome_size\tgenus\tspecies\n" > "$outfile"
            printf "%s\t%s\t%s\t%s\n" "$tax_id" "$expected_genome_size" "$genus" "$species" >> "$outfile"

            echo "$(pwd)/$outfile" >> tsv_manifest.txt
        done
    >>>

    output {
        Array[File] taxid_and_genome_size_tsvs = read_lines("tsv_manifest.txt")
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks:  "local-disk ~{disk_size_gb} SSD"
        cpu:    cpu
    }
}

task hifi_read_stats {
    meta {
        description: "Compile summary metrics from individual pipeline outputs into a single TSV, TXT, and formatted report."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
    }

    parameter_meta {
        sample_name: {
            description: "Sample name used as output prefix (e.g. bc2019)"
        }
        bam_stats_tsv: {
            description: "Output from bam_to_fastq_and_stats task",
            patterns: ["*.bam_stats.tsv"]
        }
        seqkit_stats_tsv: {
            description: "Output from hifi_seqkit_stats task (stats)",
            patterns: ["*.seqkit_stats.tsv"]
        }
        seqkit_fx2tab_tsv: {
            description: "Output from hifi_seqkit_stats task (fx2tab)",
            patterns: ["*.seqkit_fx2tab.tsv"]
        }
        kraken2_stats_tsv: {
            description: "Output from hifi_kraken2 task",
            patterns: ["*.kraken2_stats.tsv"]
        }
        taxid_and_genome_size_tsv: {
            description: "Output from batch_get_taxID_and_genome_size task",
            patterns: ["*.taxid_and_genome_size.tsv"]
        }
    }

    input {
        String sample_name
        File   bam_stats_tsv
        File   seqkit_stats_tsv
        File   seqkit_fx2tab_tsv
        File   kraken2_stats_tsv
        File   taxid_and_genome_size_tsv

        String  docker_image = "ubuntu:22.04"
        Int     disk_size_gb = 5
        Int     memory_gb    = 2
        Int     cpu          = 1
    }

    command <<<
        set -euo pipefail

        paste ~{seqkit_stats_tsv} ~{seqkit_fx2tab_tsv} ~{bam_stats_tsv} ~{kraken2_stats_tsv} ~{taxid_and_genome_size_tsv} | awk -v id="~{sample_name}" '
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
            # TXT report
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

            # TSV report
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

            # Formatted report
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
    >>>

    output {
        File hifi_read_stats_tsv    = "~{sample_name}.hifi_read_stats.tsv"
        File hifi_read_stats_txt    = "~{sample_name}.hifi_read_stats.txt"
        File hifi_read_stats_report = "~{sample_name}.hifi_read_stats.report.txt"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks:  "local-disk ~{disk_size_gb} SSD"
        cpu:    cpu
    }
}

task merge_sample_stats {
    meta {
        description: "Merge per-sample hifi_read_stats TSVs into a single multi-row TSV, adding bam and strain columns from the original sample sheet."
    }

    input {
        Array[File]   sample_tsvs
        Array[String] bam_paths
        Array[String] strains
    }

    File bam_paths_file = write_lines(bam_paths)
    File strains_file   = write_lines(strains)

    command <<<
        set -euo pipefail

        # Build merged TSV body (header from first file, data rows from all)
        head -n 1 ~{sample_tsvs[0]} > merged_body.tsv
        for f in ~{sep=' ' sample_tsvs}; do
            tail -n +2 "$f" >> merged_body.tsv
        done

        # Build extra columns file (bam, strain)
        printf "bam\tstrain\n" > extra_cols.tsv
        paste ~{bam_paths_file} ~{strains_file} >> extra_cols.tsv

        # Combine
        paste merged_body.tsv extra_cols.tsv > all_samples.hifi_read_stats.tsv
    >>>

    output {
        File merged_stats_tsv = "all_samples.hifi_read_stats.tsv"
    }

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        disks:  "local-disk 5 SSD"
        cpu:    1
    }
}

# ============================================================================
# Workflow
# ============================================================================

workflow hifi_pipeline {
    meta {
        description: "PacBio HiFi read QC pipeline. Processes BAM files through taxonomy lookup, BAM-to-FASTQ conversion, read quality metrics, genome size estimation, and taxonomic classification, then compiles per-sample and multi-sample summary reports."
        author: "tshea test"
        email:  "tshea@broadinstitute.org"
        allowNestedInputs: true
    }

    parameter_meta {
        input_tsv: {
            description: "Sample sheet TSV with header and columns: bam, genus, species, strain"
        }
        kraken2_db_hash: {
            description: "Kraken2 database hash.k2d file (pre-extracted)"
        }
        kraken2_db_opts: {
            description: "Kraken2 database opts.k2d file (pre-extracted)"
        }
        kraken2_db_taxo: {
            description: "Kraken2 database taxo.k2d file (pre-extracted)"
        }
        kraken2_confidence: {
            description: "Kraken2 confidence threshold (default 0.001)"
        }
    }

    input {
        File   input_tsv
        File   kraken2_db_hash    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/hash.k2d"
        File   kraken2_db_opts    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/opts.k2d"
        File   kraken2_db_taxo    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/taxo.k2d"
        Float  kraken2_confidence = 0.001
    }

    # ---- Parse sample sheet ------------------------------------------------------
    call parse_sample_sheet {
        input:
            sample_tsv = input_tsv
    }

    # ---- Batch taxonomy lookup (single VM for all samples) -----------------------
    call batch_get_taxID_and_genome_size {
        input:
            genera       = parse_sample_sheet.genera,
            species_list = parse_sample_sheet.species,
            barcodes     = parse_sample_sheet.barcodes
    }

    # ---- Per-sample processing ---------------------------------------------------
    scatter (idx in range(length(parse_sample_sheet.bam_paths))) {

        String bam_path = parse_sample_sheet.bam_paths[idx]
        String barcode  = parse_sample_sheet.barcodes[idx]

        # Convert BAM to FASTQ and compute accuracy/passes (single BAM download)
        call bam_to_fastq_and_stats {
            input:
                input_bam = bam_path
        }

        # Seqkit stats (depends on FASTQ)
        call hifi_seqkit_stats {
            input:
                input_fastq = bam_to_fastq_and_stats.fastq_gz
        }

        # Kraken2 classification (depends on FASTQ)
        call hifi_kraken2 {
            input:
                input_fastq     = bam_to_fastq_and_stats.fastq_gz,
                kraken2_db_hash = kraken2_db_hash,
                kraken2_db_opts = kraken2_db_opts,
                kraken2_db_taxo = kraken2_db_taxo,
                confidence      = kraken2_confidence
        }

        # Aggregate per-sample metrics
        call hifi_read_stats {
            input:
                sample_name               = barcode,
                bam_stats_tsv             = bam_to_fastq_and_stats.bam_stats_tsv,
                seqkit_stats_tsv          = hifi_seqkit_stats.seqkit_stats_tsv,
                seqkit_fx2tab_tsv         = hifi_seqkit_stats.seqkit_fx2tab_tsv,
                kraken2_stats_tsv         = hifi_kraken2.kraken2_stats,
                taxid_and_genome_size_tsv = batch_get_taxID_and_genome_size.taxid_and_genome_size_tsvs[idx]
        }
    }

    # ---- Merge all samples -------------------------------------------------------
    call merge_sample_stats {
        input:
            sample_tsvs = hifi_read_stats.hifi_read_stats_tsv,
            bam_paths   = parse_sample_sheet.bam_paths,
            strains     = parse_sample_sheet.strains
    }

    # ---- Pipeline outputs --------------------------------------------------------
    output {
        File         all_samples_stats_tsv   = merge_sample_stats.merged_stats_tsv
        Array[File]  per_sample_stats_tsv    = hifi_read_stats.hifi_read_stats_tsv
        Array[File]  per_sample_stats_txt    = hifi_read_stats.hifi_read_stats_txt
        Array[File]  per_sample_stats_report = hifi_read_stats.hifi_read_stats_report
        Array[File]  per_sample_fastq        = bam_to_fastq_and_stats.fastq_gz
        Array[File]  per_sample_kraken_report = hifi_kraken2.kraken_report
    }
}

