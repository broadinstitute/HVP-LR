version 1.0

import "../../structs/Structs.wdl"

# ============================================================================
# MMseqs2 — many-against-many protein/nucleotide sequence search and
# clustering. Same lab as foldseek (Steinegger); foldseek vendors a subset of
# mmseqs2 internals but does not expose the full mmseqs CLI. This file wraps
# the user-facing mmseqs subcommands as standalone WDL tasks.
#
# Conventions:
#  - Every task takes `String extra_args = ""` appended verbatim to the CLI.
#  - DB-mode tasks pass DBs in/out as tar.gz archives because mmseqs DBs are
#    multi-file artifacts and WDL 1.0 has no Directory type. Inside the
#    archive the DB files sit at the root with prefix `db` (e.g. `db.dbtype`,
#    `db.index`, ...). Downstream tasks auto-discover the prefix via the
#    `find ... -name '*.dbtype'` pattern; do not hard-code paths.
#  - `easy-*` tasks accept and return FASTA / TSV directly — no DB plumbing.
# ============================================================================

task MmseqsEasySearch {

    meta {
        description: "Run mmseqs easy-search end-to-end: builds query+target DBs from FASTA, runs prefilter+align, formats results as a BLAST m8 TSV. Use when you want a one-shot FASTA-in / TSV-out sequence-similarity search."

        tool:          "mmseqs easy-search"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            results_tsv: "mmseqs easy-search alignment TSV (BLAST m8 format by default). Columns determined by `format_output`.",
            num_hits:    "Number of rows in results_tsv (count of query-target hit pairs surviving the e-value cutoff)."
        }
    }

    parameter_meta {
        query_fasta:           "Query sequences (protein or nucleotide FASTA, optionally .gz)"
        target_fasta:          "Target sequences (protein or nucleotide FASTA, optionally .gz)"
        prefix:                "Basename for the results TSV (e.g. 'my_run' yields my_run.m8)"
        evalue_cutoff:         "Maximum e-value of reported hits (mmseqs -e). Default 0.001."
        sensitivity:           "Search sensitivity (mmseqs -s). 1.0=fast, 7.5=max. Default 5.7."
        format_output:         "Comma-separated mmseqs --format-output column spec. Default reports query,target,evalue,pident,bits,alnlen."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs easy-search invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_fasta
        File   target_fasta
        String prefix

        Float  evalue_cutoff = 0.001
        Float  sensitivity   = 5.7
        String format_output = "query,target,evalue,pident,bits,alnlen"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(query_fasta, "GB") + size(target_fasta, "GB")))

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

        mkdir -p tmp

        mmseqs easy-search \
            ~{query_fasta} \
            ~{target_fasta} \
            ~{prefix}.m8 \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            -s ~{sensitivity} \
            --format-output "~{format_output}" \
            ~{extra_args}

        wc -l < ~{prefix}.m8 | tr -d ' ' > nhits.txt
    >>>

    output {
        File results_tsv = "~{prefix}.m8"
        Int  num_hits    = read_int("nhits.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsEasyCluster {

    meta {
        description: "Run mmseqs easy-cluster end-to-end on a FASTA: builds a sequence DB, runs the cascaded clustering workflow, emits the cluster TSV (representative \\t member), representative-sequences FASTA, and a grouped all-sequences FASTA."

        tool:          "mmseqs easy-cluster"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. Nature Communications. 2018;9:2542."

        outputs: {
            cluster_tsv:     "Two-column TSV mapping cluster representative to each cluster member (one row per member, reps map to themselves).",
            rep_seq_fasta:   "FASTA of cluster representative sequences (one entry per cluster).",
            all_seqs_fasta:  "FASTA of every input sequence, grouped by cluster (cluster reps act as headers preceding their members).",
            num_clusters:    "Number of unique cluster representatives."
        }
    }

    parameter_meta {
        input_fasta:           "Sequences to cluster (protein or nucleotide FASTA, optionally .gz)"
        prefix:                "Basename for cluster outputs"
        min_seq_id:            "Mmseqs --min-seq-id: minimum sequence identity for cluster membership (0.0-1.0). Default 0.5."
        coverage:              "Mmseqs -c: minimum alignment coverage of the shorter sequence. Default 0.8."
        cov_mode:              "Mmseqs --cov-mode: coverage interpretation. 0=bidirectional, 1=target, 2=query. Default 0."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs easy-cluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   input_fasta
        String prefix

        Float min_seq_id = 0.5
        Float coverage   = 0.8
        Int   cov_mode   = 0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(8.0 * size(input_fasta, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p tmp

        mmseqs easy-cluster \
            ~{input_fasta} \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            --min-seq-id ~{min_seq_id} \
            -c ~{coverage} \
            --cov-mode ~{cov_mode} \
            ~{extra_args}

        cut -f1 ~{prefix}_cluster.tsv | sort -u | wc -l | tr -d ' ' > nclusters.txt
    >>>

    output {
        File cluster_tsv    = "~{prefix}_cluster.tsv"
        File rep_seq_fasta  = "~{prefix}_rep_seq.fasta"
        File all_seqs_fasta = "~{prefix}_all_seqs.fasta"
        Int  num_clusters   = read_int("nclusters.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsEasyLinclust {

    meta {
        description: "Run mmseqs easy-linclust: linear-time clustering for very large input sets. Less sensitive than easy-cluster but scales to billions of sequences. Use for first-pass NR reduction; switch to easy-cluster when sensitivity matters."

        tool:          "mmseqs easy-linclust"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. Nature Communications. 2018;9:2542."

        outputs: {
            cluster_tsv:     "Two-column TSV mapping cluster representative to each cluster member.",
            rep_seq_fasta:   "FASTA of cluster representative sequences.",
            all_seqs_fasta:  "FASTA of every input sequence, grouped by cluster.",
            num_clusters:    "Number of unique cluster representatives."
        }
    }

    parameter_meta {
        input_fasta:           "Sequences to cluster (protein or nucleotide FASTA, optionally .gz)"
        prefix:                "Basename for cluster outputs"
        min_seq_id:            "Mmseqs --min-seq-id: minimum sequence identity (0.0-1.0). Default 0.9 (typical NR reduction target)."
        coverage:              "Mmseqs -c: minimum alignment coverage of the shorter sequence. Default 0.8."
        cov_mode:              "Mmseqs --cov-mode: 0=bidirectional, 1=target, 2=query. Default 0."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs easy-linclust invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   input_fasta
        String prefix

        Float min_seq_id = 0.9
        Float coverage   = 0.8
        Int   cov_mode   = 0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(8.0 * size(input_fasta, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p tmp

        mmseqs easy-linclust \
            ~{input_fasta} \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            --min-seq-id ~{min_seq_id} \
            -c ~{coverage} \
            --cov-mode ~{cov_mode} \
            ~{extra_args}

        cut -f1 ~{prefix}_cluster.tsv | sort -u | wc -l | tr -d ' ' > nclusters.txt
    >>>

    output {
        File cluster_tsv    = "~{prefix}_cluster.tsv"
        File rep_seq_fasta  = "~{prefix}_rep_seq.fasta"
        File all_seqs_fasta = "~{prefix}_all_seqs.fasta"
        Int  num_clusters   = read_int("nclusters.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsEasyRbh {

    meta {
        description: "Run mmseqs easy-rbh: bidirectional best-hit search between two FASTA inputs. Produces only reciprocal best hits (a hit is reported only if the target's best hit back in the query set is the original query). Use for orthology calls between two protein sets."

        tool:          "mmseqs easy-rbh"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            results_tsv: "Reciprocal best-hit alignment TSV (BLAST m8).",
            num_hits:    "Number of RBH pairs."
        }
    }

    parameter_meta {
        query_fasta:           "Query sequences (FASTA)"
        target_fasta:          "Target sequences (FASTA)"
        prefix:                "Basename for the results TSV"
        evalue_cutoff:         "Maximum e-value (mmseqs -e). Default 0.001."
        sensitivity:           "Search sensitivity (mmseqs -s). Default 5.7."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs easy-rbh invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_fasta
        File   target_fasta
        String prefix

        Float evalue_cutoff = 0.001
        Float sensitivity   = 5.7

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(query_fasta, "GB") + size(target_fasta, "GB")))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p tmp

        mmseqs easy-rbh \
            ~{query_fasta} \
            ~{target_fasta} \
            ~{prefix}.m8 \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            -s ~{sensitivity} \
            ~{extra_args}

        wc -l < ~{prefix}.m8 | tr -d ' ' > nhits.txt
    >>>

    output {
        File results_tsv = "~{prefix}.m8"
        Int  num_hits    = read_int("nhits.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsEasyTaxonomy {

    meta {
        description: "Run mmseqs easy-taxonomy: assign each query sequence a taxonomic label by searching against a taxonomy-annotated reference DB (e.g. UniRef50 with taxonomy, NR with taxonomy). Outputs per-sequence taxonomy TSV, a Krona-ready report, and aggregate stats."

        tool:          "mmseqs easy-taxonomy"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Mirdita M, Steinegger M, Breitwieser F, Söding J, Levy Karin E. Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics. 2021;37(18):3029-3031."

        outputs: {
            lca_tsv:         "Per-sequence taxonomy assignment TSV (query, taxid, rank, lineage, name).",
            report_tsv:      "Krona-compatible taxonomy report.",
            tophit_aln:      "Top-hit alignment TSV (m8) supporting each taxonomy call."
        }
    }

    parameter_meta {
        query_fasta:           "Query sequences (protein or nucleotide FASTA, optionally .gz)"
        reference_db_tgz:      "Taxonomy-annotated mmseqs reference DB packaged as tar.gz (DB files are at the archive root with prefix `db` (e.g. `db.dbtype`)). Build with MmseqsCreateTaxDb or fetch with MmseqsDatabases."
        prefix:                "Basename for outputs"
        sensitivity:           "Mmseqs -s. Default 2.0 (taxonomy default; lower than search default because exhaustive search is rarely needed for taxonomy)."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs easy-taxonomy invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_fasta
        File   reference_db_tgz
        String prefix

        Float sensitivity = 2.0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(query_fasta, "GB")) + ceil(5.0 * size(reference_db_tgz, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p refdb tmp
        tar -xzf ~{reference_db_tgz} -C refdb

        RPRE=$(find refdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        RPRE=${RPRE%.dbtype}

        mmseqs easy-taxonomy \
            ~{query_fasta} \
            "${RPRE}" \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            -s ~{sensitivity} \
            ~{extra_args}
    >>>

    output {
        File lca_tsv    = "~{prefix}_lca.tsv"
        File report_tsv = "~{prefix}_report"
        File tophit_aln = "~{prefix}_tophit_aln"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsCreateDb {

    meta {
        description: "Convert a FASTA into an mmseqs sequence DB. The DB is a multi-file artifact; this task packages it as a single tar.gz so the WDL surface can pass it as one File between tasks. DB prefix inside archive: `db/db`."

        tool:          "mmseqs createdb"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            db_archive: "tar.gz of the mmseqs sequence DB. Untar with `tar -xzf` to get the DB files; DB prefix inside the archive is `db/db`."
        }
    }

    parameter_meta {
        input_fasta:           "Sequences to load into the DB (protein or nucleotide FASTA, optionally .gz)"
        prefix:                "Basename for the output archive (e.g. 'my_run' yields my_run.mmseqs_db.tar.gz)"
        extra_args:            "Additional command-line args appended verbatim to the mmseqs createdb invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   input_fasta
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_fasta, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p db

        mmseqs createdb \
            ~{input_fasta} \
            db/db \
            ~{extra_args}

        tar -czf ~{prefix}.mmseqs_db.tar.gz -C db .
    >>>

    output {
        File db_archive = "~{prefix}.mmseqs_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsCreateIndex {

    meta {
        description: "Build an mmseqs k-mer index for a sequence DB. The index speeds up downstream mmseqs search by precomputing the prefilter table. Use when the same target DB will be searched many times."

        tool:          "mmseqs createindex"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            db_archive: "tar.gz of the indexed mmseqs DB (DB + index files; DB files are at the archive root with prefix `db` (e.g. `db.dbtype`))."
        }
    }

    parameter_meta {
        db_archive:            "mmseqs DB tar.gz to index"
        prefix:                "Basename for the output archive"
        sensitivity:           "Mmseqs -s used at index time (must match the sensitivity used at search time). Default 5.7."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs createindex invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        Float sensitivity = 5.7

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(db_archive, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p db tmp
        tar -xzf ~{db_archive} -C db

        # Auto-discover DB prefix — packed archives use `-C db .` so the
        # `*.dbtype` file lands at the top level after extraction. Exclude
        # the headers-DB (`_h.dbtype`) sibling.
        DBPRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        DBPRE=${DBPRE%.dbtype}

        mmseqs createindex \
            "${DBPRE}" \
            tmp \
            --threads "${NUM_CPUS}" \
            -s ~{sensitivity} \
            ~{extra_args}

        tar -czf ~{prefix}.mmseqs_indexed_db.tar.gz -C db .
    >>>

    output {
        File db_archive_indexed = "~{prefix}.mmseqs_indexed_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsSearch {

    meta {
        description: "Run mmseqs search (DB-input/DB-output workflow): prefilter + align of a query DB against a target DB. Both inputs and the alignment output are tar.gz-packaged mmseqs DBs. Pair with MmseqsConvertAlis to surface the alignment as a TSV. Use when a target DB is searched many times — combine with MmseqsCreateIndex."

        tool:          "mmseqs search"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            aln_db_archive: "tar.gz of the alignment result DB (mmseqs-internal format). Convert with MmseqsConvertAlis."
        }
    }

    parameter_meta {
        query_db_archive:      "mmseqs query DB packaged as tar.gz"
        target_db_archive:     "mmseqs target DB packaged as tar.gz"
        prefix:                "Basename for the alignment-DB archive"
        evalue_cutoff:         "Mmseqs -e maximum e-value. Default 0.001."
        sensitivity:           "Mmseqs -s sensitivity. Default 5.7."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs search invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        String prefix

        Float evalue_cutoff = 0.001
        Float sensitivity   = 5.7

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB")))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p qdb tdb aln tmp
        tar -xzf ~{query_db_archive}  -C qdb
        tar -xzf ~{target_db_archive} -C tdb

        QPRE=$(find qdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find tdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}

        mmseqs search \
            "${QPRE}" \
            "${TPRE}" \
            aln/aln \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            -s ~{sensitivity} \
            ~{extra_args}

        tar -czf ~{prefix}.mmseqs_aln.tar.gz -C aln .
    >>>

    output {
        File aln_db_archive = "~{prefix}.mmseqs_aln.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsCluster {

    meta {
        description: "Run mmseqs cluster (DB-input cascade clustering). Lower-level than easy-cluster — accepts an existing mmseqs DB and emits a cluster-DB. Pair with MmseqsCreateDb upstream and result2tsv/createseqfiledb downstream for TSV/FASTA outputs."

        tool:          "mmseqs cluster"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. Nature Communications. 2018;9:2542."

        outputs: {
            cluster_db_archive: "tar.gz of the cluster-DB (mmseqs-internal format)."
        }
    }

    parameter_meta {
        db_archive:            "Input mmseqs sequence DB tar.gz"
        prefix:                "Basename for the cluster-DB archive"
        min_seq_id:            "Mmseqs --min-seq-id. Default 0.5."
        coverage:              "Mmseqs -c minimum alignment coverage. Default 0.8."
        cov_mode:              "Mmseqs --cov-mode. Default 0 (bidirectional)."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs cluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        Float min_seq_id = 0.5
        Float coverage   = 0.8
        Int   cov_mode   = 0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(10.0 * size(db_archive, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p db clu tmp
        tar -xzf ~{db_archive} -C db

        DBPRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        DBPRE=${DBPRE%.dbtype}

        mmseqs cluster \
            "${DBPRE}" \
            clu/clu \
            tmp \
            --threads "${NUM_CPUS}" \
            --min-seq-id ~{min_seq_id} \
            -c ~{coverage} \
            --cov-mode ~{cov_mode} \
            ~{extra_args}

        tar -czf ~{prefix}.mmseqs_cluster.tar.gz -C clu .
    >>>

    output {
        File cluster_db_archive = "~{prefix}.mmseqs_cluster.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsConvertAlis {

    meta {
        description: "Convert an mmseqs alignment-DB (output of MmseqsSearch) to a BLAST m8 TSV. Requires the query and target sequence DBs that were used in the search."

        tool:          "mmseqs convertalis"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            results_tsv: "Alignment TSV (BLAST m8 format).",
            num_hits:    "Number of rows in results_tsv."
        }
    }

    parameter_meta {
        query_db_archive:      "mmseqs query DB tar.gz (must match the query DB used in the upstream search)"
        target_db_archive:     "mmseqs target DB tar.gz (must match the target DB used in the upstream search)"
        aln_db_archive:        "mmseqs alignment-DB tar.gz (output of MmseqsSearch)"
        prefix:                "Basename for the results TSV"
        format_output:         "Comma-separated mmseqs --format-output column spec. Default reports query,target,evalue,pident,bits,alnlen."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs convertalis invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        File   aln_db_archive
        String prefix

        String format_output = "query,target,evalue,pident,bits,alnlen"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB") + size(aln_db_archive, "GB")))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p qdb tdb aln
        tar -xzf ~{query_db_archive}  -C qdb
        tar -xzf ~{target_db_archive} -C tdb
        tar -xzf ~{aln_db_archive}    -C aln

        QPRE=$(find qdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find tdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}
        APRE=$(find aln -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        APRE=${APRE%.dbtype}

        mmseqs convertalis \
            "${QPRE}" \
            "${TPRE}" \
            "${APRE}" \
            ~{prefix}.m8 \
            --threads "${NUM_CPUS}" \
            --format-output "~{format_output}" \
            ~{extra_args}

        wc -l < ~{prefix}.m8 | tr -d ' ' > nhits.txt
    >>>

    output {
        File results_tsv = "~{prefix}.m8"
        Int  num_hits    = read_int("nhits.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsCreateTsv {

    meta {
        description: "Convert an mmseqs cluster-DB to a two-column TSV (representative \\t member). Pair downstream of MmseqsCluster."

        tool:          "mmseqs createtsv"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017;35(11):1026-1028."

        outputs: {
            cluster_tsv:  "Two-column TSV (representative \\t member).",
            num_clusters: "Number of unique cluster representatives."
        }
    }

    parameter_meta {
        db_archive:            "mmseqs sequence DB tar.gz (the input DB that was clustered)"
        cluster_db_archive:    "mmseqs cluster-DB tar.gz (output of MmseqsCluster)"
        prefix:                "Basename for the cluster TSV"
        extra_args:            "Additional command-line args appended verbatim to the mmseqs createtsv invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        File   cluster_db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * (size(db_archive, "GB") + size(cluster_db_archive, "GB")))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p db clu
        tar -xzf ~{db_archive}         -C db
        tar -xzf ~{cluster_db_archive} -C clu

        DBPRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        DBPRE=${DBPRE%.dbtype}
        CPRE=$(find clu -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        CPRE=${CPRE%.dbtype}

        mmseqs createtsv \
            "${DBPRE}" \
            "${DBPRE}" \
            "${CPRE}" \
            ~{prefix}_cluster.tsv \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        cut -f1 ~{prefix}_cluster.tsv | sort -u | wc -l | tr -d ' ' > nclusters.txt
    >>>

    output {
        File cluster_tsv  = "~{prefix}_cluster.tsv"
        Int  num_clusters = read_int("nclusters.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsTaxonomy {

    meta {
        description: "Run mmseqs taxonomy on an mmseqs query DB against a taxonomy-annotated mmseqs reference DB. DB-mode counterpart to easy-taxonomy. Use when query DB is built upstream (e.g. with MmseqsCreateDb) and you want to keep results in mmseqs DB format for further processing."

        tool:          "mmseqs taxonomy"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Mirdita M, Steinegger M, Breitwieser F, Söding J, Levy Karin E. Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics. 2021;37(18):3029-3031."

        outputs: {
            tax_db_archive: "tar.gz of the taxonomy-DB (mmseqs-internal format). Convert to TSV with MmseqsCreateTsv-style postprocessing.",
            lca_tsv:        "Per-sequence taxonomy assignment TSV (query, taxid, rank, lineage, name)."
        }
    }

    parameter_meta {
        query_db_archive:      "mmseqs query DB tar.gz"
        reference_db_tgz:      "Taxonomy-annotated mmseqs reference DB tar.gz"
        prefix:                "Basename for outputs"
        sensitivity:           "Mmseqs -s. Default 2.0."
        extra_args:            "Additional command-line args appended verbatim to the mmseqs taxonomy invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   reference_db_tgz
        String prefix

        Float sensitivity = 2.0

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * (size(query_db_archive, "GB") + size(reference_db_tgz, "GB")))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p qdb refdb tax tmp
        tar -xzf ~{query_db_archive} -C qdb
        tar -xzf ~{reference_db_tgz} -C refdb

        QPRE=$(find qdb   -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        RPRE=$(find refdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' | head -1)
        RPRE=${RPRE%.dbtype}

        mmseqs taxonomy \
            "${QPRE}" \
            "${RPRE}" \
            tax/tax \
            tmp \
            --threads "${NUM_CPUS}" \
            -s ~{sensitivity} \
            ~{extra_args}

        mmseqs createtsv \
            "${QPRE}" \
            tax/tax \
            ~{prefix}_lca.tsv \
            --threads "${NUM_CPUS}"

        tar -czf ~{prefix}.mmseqs_tax.tar.gz -C tax .
    >>>

    output {
        File tax_db_archive = "~{prefix}.mmseqs_tax.tar.gz"
        File lca_tsv        = "~{prefix}_lca.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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

task MmseqsDatabases {

    meta {
        description: "Fetch a standard mmseqs reference DB (UniRef50, UniRef90, UniProtKB/Swiss-Prot, NR, NT, GTDB, etc.) via `mmseqs databases`. Output is packaged as a tar.gz with DB prefix `db/db` inside, matching the convention of MmseqsCreateDb."

        tool:          "mmseqs databases"
        tool_version:  "18-8cc5c"
        tool_url:      "https://github.com/soedinglab/MMseqs2"
        tool_citation: "Mirdita M, Steinegger M, Breitwieser F, Söding J, Levy Karin E. Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics. 2021;37(18):3029-3031."

        outputs: {
            db_archive: "tar.gz of the fetched DB; DB files are at the archive root with prefix `db` (e.g. `db.dbtype`)."
        }
    }

    parameter_meta {
        db_name:               "mmseqs DB name (e.g. 'UniRef50', 'UniRef90', 'NR', 'GTDB'). Full list: `mmseqs databases`."
        prefix:                "Basename for the output archive"
        extra_args:            "Additional command-line args appended verbatim to the mmseqs databases invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        String db_name
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"

        mkdir -p db tmp

        mmseqs databases \
            "~{db_name}" \
            db/db \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czf ~{prefix}.mmseqs_db.tar.gz -C db .
    >>>

    output {
        File db_archive = "~{prefix}.mmseqs_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            500,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/mmseqs2:18.0.0"
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
