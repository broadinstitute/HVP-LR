version 1.0

import "../../structs/Structs.wdl"

task FoldseekEasySearch {

    meta {
        description: "Run foldseek easy-search to find structural homologs of each query against the target set. Accepts PDB/mmCIF inputs as two arrays of files, stages them into query/ and target/ directories, and runs foldseek's full easy-search workflow (createdb on both sides, structural prefilter, alignment, format). Emits the alignment hits as a TSV plus the hit count."

        tool:          "foldseek"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            results_tsv: "Foldseek easy-search alignment hits, one row per (query, target) hit. Columns determined by `format_output`.",
            num_hits:    "Number of rows in results_tsv (count of query-target hit pairs surviving the e-value cutoff)."
        }
    }

    parameter_meta {
        query_structures:      "Per-structure files (PDB or mmCIF, optionally .gz) used as queries. Empty array is not supported."
        target_structures:     "Per-structure files (PDB or mmCIF, optionally .gz) used as the search target set. Empty array is not supported."
        prefix:                "Basename for the results TSV (e.g. 'my_run' yields my_run.m8)."
        evalue_cutoff:         "Maximum e-value of reported hits (foldseek -e). Default 0.001."
        format_output:         "Comma-separated foldseek --format-output column spec. Default reports query/target/evalue/bits/fident/alnlen/prob."
        extra_args:            "Additional command-line args appended verbatim to the foldseek easy-search invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] query_structures
        Array[File] target_structures
        String      prefix

        Float  evalue_cutoff = 0.001
        String format_output = "query,target,evalue,bits,fident,alnlen,prob"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(query_structures, "GB") + size(target_structures, "GB")))

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

        # Stage queries and targets into stable directories. Cromwell localizes
        # array files into per-input subdirectories that foldseek's directory
        # scanner doesn't recurse into reliably, so flatten via symlink.
        mkdir -p query target tmp
        for f in ~{sep=" " query_structures};  do ln -sf "$f" query/;  done
        for f in ~{sep=" " target_structures}; do ln -sf "$f" target/; done

        foldseek easy-search \
            query/ \
            target/ \
            ~{prefix}.m8 \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            --format-output "~{format_output}" \
            ~{extra_args}

        # Hit count for downstream scalar surface.
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
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekEasyCluster {

    meta {
        description: "Run foldseek easy-cluster end-to-end on a set of protein structures: builds a structure DB, runs the cascade clustering workflow, and emits the cluster TSV (representative \t member), representative sequences FASTA, and all-vs-all sequences FASTA."

        tool:          "foldseek easy-cluster"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            cluster_tsv:     "Two-column TSV mapping cluster representative to each cluster member (one row per member, reps map to themselves).",
            rep_seq_fasta:   "FASTA of cluster representative sequences (one entry per cluster).",
            all_seqs_fasta:  "FASTA of every input sequence, grouped by cluster (cluster reps act as headers preceding their members).",
            num_clusters:    "Number of unique cluster representatives."
        }
    }

    parameter_meta {
        structures:            "Per-structure files (PDB or mmCIF, optionally .gz) to cluster"
        prefix:                "Basename for cluster outputs"
        min_seq_id:            "Foldseek --min-seq-id: minimum sequence identity for cluster membership. Default 0.0 (structure-only)."
        coverage:              "Foldseek -c: minimum alignment coverage of the shorter sequence. Default 0.8."
        extra_args:            "Additional command-line args appended verbatim to the foldseek easy-cluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] structures
        String      prefix

        Float min_seq_id = 0.0
        Float coverage   = 0.8

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(structures, "GB"))

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

        mkdir -p input tmp
        for f in ~{sep=" " structures}; do ln -sf "$f" input/; done

        foldseek easy-cluster \
            input/ \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            --min-seq-id ~{min_seq_id} \
            -c ~{coverage} \
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
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekEasyRbh {

    meta {
        description: "Run foldseek easy-rbh to find reciprocal best hits between a query and a target set of protein structures end-to-end (createdb on both sides, structural prefilter, alignment, RBH selection, format)."

        tool:          "foldseek easy-rbh"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            rbh_tsv:  "Reciprocal best hits, one row per (query, target) RBH pair.",
            num_hits: "Number of rows in rbh_tsv."
        }
    }

    parameter_meta {
        query_structures:      "Per-structure query files (PDB or mmCIF, optionally .gz)"
        target_structures:     "Per-structure target files (PDB or mmCIF, optionally .gz)"
        prefix:                "Basename for the RBH TSV"
        evalue_cutoff:         "Foldseek -e maximum e-value for considered hits. Default 0.001."
        format_output:         "Foldseek --format-output column spec. Default reports query/target/evalue/bits/fident/alnlen/prob."
        extra_args:            "Additional command-line args appended verbatim to the foldseek easy-rbh invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] query_structures
        Array[File] target_structures
        String      prefix

        Float  evalue_cutoff = 0.001
        String format_output = "query,target,evalue,bits,fident,alnlen,prob"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(query_structures, "GB") + size(target_structures, "GB")))

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

        mkdir -p query target tmp
        for f in ~{sep=" " query_structures};  do ln -sf "$f" query/;  done
        for f in ~{sep=" " target_structures}; do ln -sf "$f" target/; done

        foldseek easy-rbh \
            query/ \
            target/ \
            ~{prefix}.rbh.m8 \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            --format-output "~{format_output}" \
            ~{extra_args}

        wc -l < ~{prefix}.rbh.m8 | tr -d ' ' > nhits.txt
    >>>

    output {
        File rbh_tsv  = "~{prefix}.rbh.m8"
        Int  num_hits = read_int("nhits.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekEasyMultimerSearch {

    meta {
        description: "Run foldseek easy-multimersearch to search query multimer structures against a set of target multimer structures end-to-end. Surfaces both the per-chain alignment TSV and the per-multimer report."

        tool:          "foldseek easy-multimersearch"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            results_tsv: "Per-chain alignment hits.",
            report_tsv:  "Per-multimer report aggregating chain-level hits into complex-level matches."
        }
    }

    parameter_meta {
        query_structures:      "Per-structure query files of multimer complexes (PDB or mmCIF, optionally .gz)"
        target_structures:     "Per-structure target files of multimer complexes (PDB or mmCIF, optionally .gz)"
        prefix:                "Basename for output files"
        extra_args:            "Additional command-line args appended verbatim to the foldseek easy-multimersearch invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] query_structures
        Array[File] target_structures
        String      prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(query_structures, "GB") + size(target_structures, "GB")))

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

        mkdir -p query target tmp
        for f in ~{sep=" " query_structures};  do ln -sf "$f" query/;  done
        for f in ~{sep=" " target_structures}; do ln -sf "$f" target/; done

        foldseek easy-multimersearch \
            query/ \
            target/ \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}
    >>>

    output {
        File results_tsv = "~{prefix}"
        File report_tsv  = "~{prefix}_report"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekEasyMultimerCluster {

    meta {
        description: "Run foldseek easy-multimercluster end-to-end on a set of multimer structures: builds the structure DB, runs multimer-level cascade clustering, and emits the cluster TSV plus a per-complex info TSV."

        tool:          "foldseek easy-multimercluster"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            cluster_tsv:    "Cluster representative -> member TSV at multimer (complex) level.",
            rep_seq_fasta:  "FASTA of cluster representative sequences (one entry per cluster)."
        }
    }

    parameter_meta {
        structures:            "Per-structure files of multimer complexes (PDB or mmCIF, optionally .gz)"
        prefix:                "Basename for output files"
        extra_args:            "Additional command-line args appended verbatim to the foldseek easy-multimercluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] structures
        String      prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(structures, "GB"))

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

        mkdir -p input tmp
        for f in ~{sep=" " structures}; do ln -sf "$f" input/; done

        foldseek easy-multimercluster \
            input/ \
            ~{prefix} \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}
    >>>

    output {
        File cluster_tsv    = "~{prefix}_cluster.tsv"
        File rep_seq_fasta  = "~{prefix}_rep_seq.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCreateDb {

    meta {
        description: "Convert PDB/mmCIF (or .gz) structure files into a foldseek structure DB. The DB is a multi-file artifact; this task packages it as a single tar.gz so the WDL surface can pass it as one File between tasks."

        tool:          "foldseek createdb"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            db_archive: "tar.gz of the foldseek structure DB. Untar with `tar -xzf` to get the DB files; DB prefix inside the archive is `db/db`."
        }
    }

    parameter_meta {
        structures:            "Per-structure files (PDB or mmCIF, optionally .gz) to load into the DB"
        prefix:                "Basename for the output archive (e.g. 'my_run' yields my_run.foldseek_db.tar.gz)"
        extra_args:            "Additional command-line args appended verbatim to the foldseek createdb invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File] structures
        String      prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(structures, "GB"))

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

        mkdir -p input db
        for f in ~{sep=" " structures}; do ln -sf "$f" input/; done

        foldseek createdb \
            input/ \
            db/db \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_db.tar.gz -C db .
    >>>

    output {
        File db_archive = "~{prefix}.foldseek_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCreateDbFromFasta {

    meta {
        description: "Build a foldseek structure DB directly from an amino-acid FASTA by predicting 3Di structural tokens with ProstT5 (foldseek's bundled AA->3Di language model). Output is a tar.gz-packaged foldseek DB that downstream FoldseekSearch / FoldseekEasySearch consume as either side. Avoids needing an external structure predictor (ESMFold etc.) when the goal is structural search of a protein set against an existing structure DB like BFVD."

        tool:          "foldseek createdb"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "Heinzinger M, Weissenow K, Sanchez JG, et al. Bilingual language model for protein sequence and structure. NAR Genomics and Bioinformatics. 2024;6(4):lqae150."

        outputs: {
            db_archive: "tar.gz of the foldseek structure DB built from AA via ProstT5. Untar to get the DB files at the archive root (prefix `db`)."
        }
    }

    parameter_meta {
        protein_fasta:         "Amino-acid FASTA (optionally .gz) to convert into a foldseek structure DB via ProstT5"
        prefix:                "Basename for the output archive (e.g. 'my_run' yields my_run.foldseek_prostt5_db.tar.gz)"
        prostt5_weights_tgz:   "ProstT5 model weights packaged as tar.gz. Source: https://foldseek.steineggerlab.workers.dev/prostt5-f16-gguf.tar.gz (~2.1 GB; single `prostt5-f16.gguf` at archive root). Stage once in Terra workspace data; reuse across runs."
        use_gpu:               "When true (default), run ProstT5 inference on a CUDA device using the foldseek-gpu image with --gpu 1 and a single nvidia-tesla-t4. When false, run on CPU using the foldseek image with no GPU attached. GPU is ~10-100x faster on large protein sets."
        extra_args:            "Additional command-line args appended verbatim to the foldseek createdb invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   protein_fasta
        String prefix
        File   prostt5_weights_tgz

        Boolean use_gpu = true

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(15.0 * size(protein_fasta, "GB")) + ceil(3.0 * size(prostt5_weights_tgz, "GB"))
    String default_docker = if use_gpu
        then "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.1"
        else "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"

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

        # ---- GPU diagnostics (use_gpu=true only) ----
        # Surfaces driver/runtime state so a CPU-speed run on a GPU-declared
        # task is debuggable from the log alone. nvcc is absent from the
        # foldseek-gpu runtime image (toolkit only in the builder stage) —
        # NO_NVCC is expected and not a failure.
        if ~{if use_gpu then "true" else "false"}; then
            echo "=== GPU diag ==="
            nvidia-smi || echo "NO_NVIDIA_SMI (driver missing or no GPU attached)"
            nvcc --version 2>/dev/null || echo "NO_NVCC (expected on runtime image)"
            echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
            echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-unset}"
            ls -1 /usr/local/cuda*/lib64/libcudart* 2>/dev/null || echo "NO_LIBCUDART_FOUND"
            echo "=== end GPU diag ==="
        fi

        # Stage ProstT5 weights. Tarball has a single `prostt5-f16.gguf` at
        # the archive root. Foldseek's --prostt5-model takes the .gguf FILE
        # path for GGUF mode (vs the directory path for the SafeTensors HF
        # bundle, which the upstream README documents but we don't use).
        # holding the .gguf, not the file itself.
        mkdir -p prostt5 db tmp
        tar -xzf ~{prostt5_weights_tgz} -C prostt5

        # Decompress AA FASTA if gzipped — foldseek createdb accepts plain
        # FASTA for the ProstT5 path. Uniform local name simplifies the CLI.
        if [[ "~{protein_fasta}" == *.gz ]]; then
            gunzip -c ~{protein_fasta} > input.faa
        else
            cp ~{protein_fasta} input.faa
        fi

        foldseek createdb \
            input.faa \
            db/db \
            --prostt5-model prostt5/prostt5-f16.gguf \
            --threads "${NUM_CPUS}" \
            ~{if use_gpu then "--gpu 1" else ""} \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_prostt5_db.tar.gz -C db .
    >>>

    output {
        File db_archive = "~{prefix}.foldseek_prostt5_db.tar.gz"
    }

    #########################
    # When use_gpu=true: ProstT5 inference on T4 (10-100x faster than CPU).
    # gpuType is a single value (Cromwell GCP backend has no multi-type fallback).
    # All zones MUST be in a single region: GCP Batch (Cromwell's backend on
    # Terra) rejects multi-region zone lists with "allocation_policy field is
    # invalid. all specified locations end up in more than one regions".
    # us-central1-c/-f keep T4 capacity within the workspace's home region.
    # When use_gpu=false: gpuCount=0 disables GPU attach; zones still
    # constrain to us-central1 to stay collocated with the GPU path's region.
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             default_docker
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
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               if use_gpu then 1 else 0
        zones:                  "us-central1-c us-central1-f"
    }
}

task FoldseekSearch {

    meta {
        description: "Run foldseek search (the DB-input/DB-output workflow): structural prefilter + alignment of a query DB against a target DB. Both inputs and the alignment output are tar.gz-packaged foldseek DBs. Pair with FoldseekConvertAlis to surface the alignment as a TSV."

        tool:          "foldseek search"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            aln_db_archive: "tar.gz of the alignment result DB (foldseek-internal format). Convert with FoldseekConvertAlis."
        }
    }

    parameter_meta {
        query_db_archive:      "Foldseek query DB packaged as tar.gz (e.g. output of FoldseekCreateDb). When use_gpu=true, query was built with FoldseekCreateDbFromFasta use_gpu=true."
        target_db_archive:     "Foldseek target DB packaged as tar.gz. When use_gpu=true, MUST be a PADDED foldseek DB (see makepaddedseqdb recipe in inline command comment block)."
        prefix:                "Basename for the alignment-DB archive"
        evalue_cutoff:         "Foldseek -e maximum e-value. Default 0.001."
        use_gpu:               "When true (default), engage CUDA prefilter with --gpu 1 --prefilter-mode 1 on the foldseek-gpu image with a single nvidia-tesla-t4. Requires target_db_archive to be a padded DB. When false, run CPU-only search on the foldseek image with no GPU attached."
        extra_args:            "Additional command-line args appended verbatim to the foldseek search invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        String prefix

        Float evalue_cutoff = 0.001

        Boolean use_gpu = true

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB")))
    String default_docker = if use_gpu
        then "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.1"
        else "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"

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

        # ---- GPU diagnostics (use_gpu=true only) ----
        # Surfaces driver/runtime state so a CPU-speed run on a GPU-declared
        # task is debuggable from the log alone. nvcc is absent from the
        # foldseek-gpu runtime image (toolkit only in the builder stage) —
        # NO_NVCC is expected and not a failure.
        if ~{if use_gpu then "true" else "false"}; then
            echo "=== GPU diag ==="
            nvidia-smi || echo "NO_NVIDIA_SMI (driver missing or no GPU attached)"
            nvcc --version 2>/dev/null || echo "NO_NVCC (expected on runtime image)"
            echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
            echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-unset}"
            ls -1 /usr/local/cuda*/lib64/libcudart* 2>/dev/null || echo "NO_LIBCUDART_FOUND"
            echo "=== end GPU diag ==="
        fi

        mkdir -p querydb targetdb result tmp
        tar -xzf ~{query_db_archive}  -C querydb
        tar -xzf ~{target_db_archive} -C targetdb

        # Query side is never padded (createdb has no padding mode) — pick
        # the unpadded primary, explicitly skipping any _pad* siblings that
        # may have been bundled in alongside.
        QPRE=$(find querydb  -maxdepth 2 -name '*.dbtype' \
                 ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' \
                 ! -name '*_pad*.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}

        # Target side: when use_gpu=true REQUIRE a *_pad.dbtype prefix
        # (built via makepaddedseqdb — see comment block below). When false,
        # use the unpadded primary and explicitly skip any _pad* sibling
        # so a padded tarball still works in CPU mode.
        if ~{if use_gpu then "true" else "false"}; then
            TPRE=$(find targetdb -maxdepth 2 -name '*_pad.dbtype' | head -1)
            if [[ -z "${TPRE}" ]]; then
                echo "ERROR: use_gpu=true requires a padded foldseek DB; no *_pad.dbtype found in target_db_archive." >&2
                echo "       See the BFVD padding recipe in the FoldseekSearch task comment block." >&2
                exit 1
            fi
        else
            TPRE=$(find targetdb -maxdepth 2 -name '*.dbtype' \
                     ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' \
                     ! -name '*_pad*.dbtype' | head -1)
        fi
        TPRE=${TPRE%.dbtype}
        echo "Selected query prefix:  ${QPRE}"
        echo "Selected target prefix: ${TPRE}"

        # -a stores the backtrace cigar in the alignment DB so downstream
        # FoldseekConvertAlis can format the hits without recomputing.
        # Override with --alignment-mode 0 via extra_args to disable.
        #
        # When use_gpu=true: --gpu 1 --prefilter-mode 1 engages CUDA prefilter
        # (foldseek-gpu image + nvidia-tesla-t4 default in runtime block).
        # When use_gpu=false: CPU-only search; flags omitted; foldseek image.
        #
        # IMPORTANT (use_gpu=true only): target_db_archive MUST contain a
        # PADDED foldseek DB at the `<prefix>_pad` sibling of the primary
        # `<prefix>`. The shipping BFVD tarball
        # (gs://hvp-tech-dev/databases/foldseek/bfvd_foldseekdb.tar.gz) is
        # unpadded — passing it to GPU search will fail at TPRE detection
        # below. Build a self-contained padded tarball once on a GPU
        # workstation:
        #
        # Foldseek 10's `makepaddedseqdb` is self-contained: it pads the 3Di
        # (_ss) DB and emits `<prefix>_pad`, `<prefix>_pad_h`, `<prefix>_pad_ca`
        # as SOFTLINKS to the originals via `renamedbkeys --subdb-mode 1`.
        # Those softlinks use the container's absolute mount path (e.g.
        # /w/db/bfvd_h), so tarring OUTSIDE the container captures dangling
        # links — the resulting tarball is tiny (only the real `_pad_ss`
        # file survives) and structurally unusable. Run BOTH the pad and the
        # tar INSIDE the container so `tar -h` dereferences the softlinks
        # to real bytes:
        #
        #     gcloud storage cp \
        #         gs://hvp-tech-dev/databases/foldseek/bfvd_foldseekdb.tar.gz .
        #     mkdir db && tar -xzf bfvd_foldseekdb.tar.gz -C db
        #     docker run --rm --gpus all -v "$PWD:/w" -w /w/db \
        #         --entrypoint bash \
        #         us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek-gpu:10.0.1 \
        #         -c '
        #             set -euxo pipefail
        #             PRE=$(find . -maxdepth 1 -name "*.dbtype" \
        #                     ! -name "*_h.dbtype" ! -name "*_ss.dbtype" ! -name "*_ca.dbtype" | head -1)
        #             PRE=${PRE%.dbtype}
        #             /usr/local/bin/foldseek makepaddedseqdb "${PRE}" "${PRE}_pad"
        #             # -h dereferences the _pad / _pad_h / _pad_ca softlinks
        #             tar -czhf /w/bfvd_foldseekdb_padded.tar.gz -C /w/db .
        #         '
        #     gcloud storage cp bfvd_foldseekdb_padded.tar.gz \
        #         gs://hvp-tech-dev/databases/foldseek/
        #
        # Then point the workspace var `bfvd_db_tgz` at the `_padded.tar.gz`.
        # Resulting tarball is roughly the SAME size as the unpadded one
        # (~14 GB for BFVD) — if dramatically smaller, the softlinks were
        # not dereferenced and the tarball is unusable. When use_gpu=false
        # the unpadded BFVD tarball works directly.
        foldseek search \
            "${QPRE}" \
            "${TPRE}" \
            result/alnDB \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            -a \
            ~{if use_gpu then "--gpu 1 --prefilter-mode 1" else ""} \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_aln_db.tar.gz -C result .
    >>>

    output {
        File aln_db_archive = "~{prefix}.foldseek_aln_db.tar.gz"
    }

    #########################
    # When use_gpu=true: foldseek search --gpu prefilter dramatically faster
    # vs large structure DBs (BFVD ~350k targets). Single gpuType per task;
    # `zones` constrained to a single region (us-central1) per GCP Batch's
    # allocation_policy rule (Cromwell rejects multi-region zone lists).
    # When use_gpu=false: gpuCount=0 detaches the accelerator; zones still
    # held to us-central1 to stay collocated with the GPU path's region.
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             default_docker
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
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               if use_gpu then 1 else 0
        zones:                  "us-central1-c us-central1-f"
    }
}

task FoldseekRbh {

    meta {
        description: "Run foldseek rbh (DB-input/DB-output reciprocal best hits) between a query DB and a target DB. Pair with FoldseekConvertAlis to surface the result as a TSV."

        tool:          "foldseek rbh"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            rbh_db_archive: "tar.gz of the RBH alignment DB. Convert with FoldseekConvertAlis."
        }
    }

    parameter_meta {
        query_db_archive:      "Foldseek query DB packaged as tar.gz"
        target_db_archive:     "Foldseek target DB packaged as tar.gz"
        prefix:                "Basename for the RBH-DB archive"
        evalue_cutoff:         "Foldseek -e maximum e-value. Default 0.001."
        extra_args:            "Additional command-line args appended verbatim to the foldseek rbh invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        String prefix

        Float evalue_cutoff = 0.001

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB")))

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

        mkdir -p querydb targetdb result tmp
        tar -xzf ~{query_db_archive}  -C querydb
        tar -xzf ~{target_db_archive} -C targetdb

        QPRE=$(find querydb  -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find targetdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}

        # -a stores backtrace cigar so downstream FoldseekConvertAlis works.
        foldseek rbh \
            "${QPRE}" \
            "${TPRE}" \
            result/rbhDB \
            tmp \
            --threads "${NUM_CPUS}" \
            -e ~{evalue_cutoff} \
            -a \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_rbh_db.tar.gz -C result .
    >>>

    output {
        File rbh_db_archive = "~{prefix}.foldseek_rbh_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCluster {

    meta {
        description: "Run foldseek cluster (DB-input/DB-output cascade clustering) on a structure DB. Output is a cluster DB; convert with createtsv or pair with downstream tasks."

        tool:          "foldseek cluster"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            cluster_db_archive: "tar.gz of the foldseek cluster DB."
        }
    }

    parameter_meta {
        db_archive:            "Foldseek structure DB packaged as tar.gz"
        prefix:                "Basename for the cluster-DB archive"
        min_seq_id:            "Foldseek --min-seq-id: minimum sequence identity. Default 0.0 (structure-only)."
        coverage:              "Foldseek -c: minimum coverage of the shorter sequence. Default 0.8."
        extra_args:            "Additional command-line args appended verbatim to the foldseek cluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        Float min_seq_id = 0.0
        Float coverage   = 0.8

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(db_archive, "GB"))

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

        mkdir -p db result tmp
        tar -xzf ~{db_archive} -C db

        PRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        PRE=${PRE%.dbtype}

        foldseek cluster \
            "${PRE}" \
            result/clusterDB \
            tmp \
            --threads "${NUM_CPUS}" \
            --min-seq-id ~{min_seq_id} \
            -c ~{coverage} \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_cluster_db.tar.gz -C result .
    >>>

    output {
        File cluster_db_archive = "~{prefix}.foldseek_cluster_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekMultimerSearch {

    meta {
        description: "Run foldseek multimersearch (DB-input/DB-output multimer-level search). Both inputs and the output are tar.gz-packaged foldseek DBs."

        tool:          "foldseek multimersearch"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            multimer_db_archive: "tar.gz of the multimer alignment DB."
        }
    }

    parameter_meta {
        query_db_archive:      "Foldseek query multimer DB packaged as tar.gz"
        target_db_archive:     "Foldseek target multimer DB packaged as tar.gz"
        prefix:                "Basename for the multimer-DB archive"
        extra_args:            "Additional command-line args appended verbatim to the foldseek multimersearch invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB")))

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

        mkdir -p querydb targetdb result tmp
        tar -xzf ~{query_db_archive}  -C querydb
        tar -xzf ~{target_db_archive} -C targetdb

        QPRE=$(find querydb  -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find targetdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}

        # -a stores backtrace cigar so downstream conversion tasks work.
        foldseek multimersearch \
            "${QPRE}" \
            "${TPRE}" \
            result/multimerDB \
            tmp \
            --threads "${NUM_CPUS}" \
            -a \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_multimer_db.tar.gz -C result .
    >>>

    output {
        File multimer_db_archive = "~{prefix}.foldseek_multimer_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekMultimerCluster {

    meta {
        description: "Run foldseek multimercluster (DB-input/DB-output multimer-level clustering) on a multimer DB."

        tool:          "foldseek multimercluster"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            multimer_cluster_db_archive: "tar.gz of the multimer cluster DB."
        }
    }

    parameter_meta {
        db_archive:            "Foldseek multimer DB packaged as tar.gz"
        prefix:                "Basename for the multimer-cluster-DB archive"
        extra_args:            "Additional command-line args appended verbatim to the foldseek multimercluster invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(db_archive, "GB"))

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

        mkdir -p db result tmp
        tar -xzf ~{db_archive} -C db

        PRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        PRE=${PRE%.dbtype}

        foldseek multimercluster \
            "${PRE}" \
            result/multimerClusterDB \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_multimer_cluster_db.tar.gz -C result .
    >>>

    output {
        File multimer_cluster_db_archive = "~{prefix}.foldseek_multimer_cluster_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekDatabases {

    meta {
        description: "Download a foldseek-hosted reference DB by name (e.g. PDB, AFDB50, AFDB-SwissProt, ProstT5) into a tar.gz-packaged structure DB usable by downstream foldseek tasks. Uses the foldseek-bundled aria2 downloader."

        tool:          "foldseek databases"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            db_archive: "tar.gz of the downloaded foldseek DB."
        }
    }

    parameter_meta {
        db_name:               "Foldseek-hosted DB name. Run `foldseek databases` with no args inside the image to see the current list."
        prefix:                "Basename for the output archive"
        extra_args:            "Additional command-line args appended verbatim to the foldseek databases invocation"
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

        mkdir -p db tmp

        foldseek databases \
            ~{db_name} \
            db/db \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_db.tar.gz -C db .
    >>>

    output {
        File db_archive = "~{prefix}.foldseek_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            500,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCreateIndex {

    meta {
        description: "Precompute a foldseek index for a structure DB to speed up downstream searches against it. Output is the original DB plus the index, repackaged into a single tar.gz."

        tool:          "foldseek createindex"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            indexed_db_archive: "tar.gz of the DB plus its precomputed foldseek index."
        }
    }

    parameter_meta {
        db_archive:            "Foldseek structure DB packaged as tar.gz"
        prefix:                "Basename for the indexed-DB archive"
        extra_args:            "Additional command-line args appended verbatim to the foldseek createindex invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(db_archive, "GB"))

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

        mkdir -p db tmp
        tar -xzf ~{db_archive} -C db

        PRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        PRE=${PRE%.dbtype}

        foldseek createindex \
            "${PRE}" \
            tmp \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czf ~{prefix}.foldseek_indexed_db.tar.gz -C db .
    >>>

    output {
        File indexed_db_archive = "~{prefix}.foldseek_indexed_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCreateCluSearchDb {

    meta {
        description: "Build a searchable cluster DB from an existing structure DB and a foldseek cluster DB, enabling faster downstream searches that respect cluster representatives."

        tool:          "foldseek createclusearchdb"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            cluster_search_db_archive: "tar.gz of the searchable cluster DB."
        }
    }

    parameter_meta {
        db_archive:            "Foldseek structure DB packaged as tar.gz"
        cluster_db_archive:    "Foldseek cluster DB packaged as tar.gz (output of FoldseekCluster)"
        prefix:                "Basename for the cluster-search-DB archive"
        extra_args:            "Additional command-line args appended verbatim to the foldseek createclusearchdb invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        File   cluster_db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(db_archive, "GB") + size(cluster_db_archive, "GB")))

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

        mkdir -p db clu result tmp
        tar -xzf ~{db_archive}         -C db
        tar -xzf ~{cluster_db_archive} -C clu

        PRE=$(find db  -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        PRE=${PRE%.dbtype}
        CPRE=$(find clu -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        CPRE=${CPRE%.dbtype}

        # createclusearchdb creates a symlink-DB pointing back at the source
        # structure DB; the symlink targets break unless both prefixes are
        # absolute paths. Resolve before invoking.
        PRE_ABS=$(realpath "${PRE}.dbtype"); PRE_ABS=${PRE_ABS%.dbtype}
        CPRE_ABS=$(realpath "${CPRE}.dbtype"); CPRE_ABS=${CPRE_ABS%.dbtype}
        RESULT_ABS=$(realpath result)/cluSearchDB

        foldseek createclusearchdb \
            "${PRE_ABS}" \
            "${CPRE_ABS}" \
            "${RESULT_ABS}" \
            --threads "${NUM_CPUS}" \
            ~{extra_args}

        tar -czhf ~{prefix}.foldseek_clu_search_db.tar.gz -C result .
    >>>

    output {
        File cluster_search_db_archive = "~{prefix}.foldseek_clu_search_db.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekConvertAlis {

    meta {
        description: "Convert a foldseek alignment DB (output of foldseek search / rbh / multimersearch) into a BLAST-tab, SAM, or custom-format TSV consumable by downstream tools."

        tool:          "foldseek convertalis"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            results_tsv: "Formatted alignment results, one row per hit. Columns determined by `format_output`.",
            num_hits:    "Number of rows in results_tsv."
        }
    }

    parameter_meta {
        query_db_archive:      "Foldseek query DB packaged as tar.gz (same DB used to produce aln_db_archive)"
        target_db_archive:     "Foldseek target DB packaged as tar.gz (same DB used to produce aln_db_archive)"
        aln_db_archive:        "Foldseek alignment DB packaged as tar.gz (output of FoldseekSearch / FoldseekRbh / FoldseekMultimerSearch)"
        prefix:                "Basename for the output TSV"
        format_output:         "Foldseek --format-output column spec. Default reports query/target/evalue/bits/fident/alnlen/prob."
        extra_args:            "Additional command-line args appended verbatim to the foldseek convertalis invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        File   aln_db_archive
        String prefix

        String format_output = "query,target,evalue,bits,fident,alnlen,prob"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB") + size(aln_db_archive, "GB")))

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

        mkdir -p querydb targetdb alndb
        tar -xzf ~{query_db_archive}  -C querydb
        tar -xzf ~{target_db_archive} -C targetdb
        tar -xzf ~{aln_db_archive}    -C alndb

        QPRE=$(find querydb  -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find targetdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}
        APRE=$(find alndb    -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        APRE=${APRE%.dbtype}

        foldseek convertalis \
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
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekConvert2Pdb {

    meta {
        description: "Convert a foldseek structure DB back into PDB files. Default mode produces a single multi-model PDB; set --multi-model 0 via extra_args to emit one PDB per entry (tarred into pdb_output)."

        tool:          "foldseek convert2pdb"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            pdb_output: "Either a single multi-model PDB file (default) or a tar.gz of per-entry PDBs, depending on foldseek's --multi-model flag passed through extra_args."
        }
    }

    parameter_meta {
        db_archive:            "Foldseek structure DB packaged as tar.gz"
        prefix:                "Basename for the output"
        extra_args:            "Additional command-line args appended verbatim to the foldseek convert2pdb invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * size(db_archive, "GB"))

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

        mkdir -p db pdbout
        tar -xzf ~{db_archive} -C db

        PRE=$(find db -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        PRE=${PRE%.dbtype}

        foldseek convert2pdb \
            "${PRE}" \
            pdbout/~{prefix} \
            ~{extra_args}

        # Foldseek writes either a single .pdb file at pdbout/<prefix> or a
        # directory of per-entry PDBs at pdbout/<prefix>/. Surface both shapes
        # under a single deterministic output name.
        if [[ -d pdbout/~{prefix} ]]; then
            tar -czf ~{prefix}.pdb.tar.gz -C pdbout/~{prefix} .
            mv ~{prefix}.pdb.tar.gz ~{prefix}.pdb_output
        else
            mv pdbout/~{prefix} ~{prefix}.pdb_output
        fi
    >>>

    output {
        File pdb_output = "~{prefix}.pdb_output"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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

task FoldseekCreateMultimerReport {

    meta {
        description: "Convert a foldseek multimer alignment DB into a TSV report aggregating chain-level hits into complex-level matches."

        tool:          "foldseek createmultimerreport"
        tool_version:  "10.0.1"
        tool_url:      "https://github.com/steineggerlab/foldseek"
        tool_citation: "van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology. 2024;42(2):243-246."

        outputs: {
            report_tsv: "Per-multimer alignment report TSV."
        }
    }

    parameter_meta {
        query_db_archive:      "Foldseek query multimer DB packaged as tar.gz"
        target_db_archive:     "Foldseek target multimer DB packaged as tar.gz"
        complex_db_archive:    "Foldseek multimer alignment DB packaged as tar.gz (output of FoldseekMultimerSearch)"
        prefix:                "Basename for the output TSV"
        extra_args:            "Additional command-line args appended verbatim to the foldseek createmultimerreport invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   query_db_archive
        File   target_db_archive
        File   complex_db_archive
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(10.0 * (size(query_db_archive, "GB") + size(target_db_archive, "GB") + size(complex_db_archive, "GB")))

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

        mkdir -p querydb targetdb complexdb
        tar -xzf ~{query_db_archive}   -C querydb
        tar -xzf ~{target_db_archive}  -C targetdb
        tar -xzf ~{complex_db_archive} -C complexdb

        QPRE=$(find querydb   -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        QPRE=${QPRE%.dbtype}
        TPRE=$(find targetdb  -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        TPRE=${TPRE%.dbtype}
        CPRE=$(find complexdb -maxdepth 2 -name '*.dbtype' ! -name '*_h.dbtype' ! -name '*_ca.dbtype' ! -name '*_ss.dbtype' | head -1)
        CPRE=${CPRE%.dbtype}

        foldseek createmultimerreport \
            "${QPRE}" \
            "${TPRE}" \
            "${CPRE}" \
            ~{prefix}.multimer_report.tsv \
            --threads "${NUM_CPUS}" \
            ~{extra_args}
    >>>

    output {
        File report_tsv = "~{prefix}.multimer_report.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/foldseek:10.0.1"
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
