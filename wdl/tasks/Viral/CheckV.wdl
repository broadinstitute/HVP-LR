version 1.0

import "../../structs/Structs.wdl"

task CheckV {

    meta {
        description: "Run CheckV end-to-end on a viral FASTA to assess completeness and contamination. Spaces in FASTA headers are replaced with underscores before running to prevent CheckV from producing duplicate IDs (it truncates headers at whitespace, causing downstream collisions in the summary script). Called once per upstream viral detection tool (geNomad, VirSorter2)."

        tool:         "CheckV"
        tool_version: "1.0.3"
        tool_url:     "https://bitbucket.org/berkeleylab/checkv"
        tool_citation: "Nayfach S, et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol. 2021;39(5):578-585."

        outputs: {
            quality_summary_tsv: "CheckV quality_summary.tsv: per-sequence completeness, contamination, provirus status, and MIUVIG quality tier",
            viruses_fna:         "CheckV viruses.fna: sequences classified as complete or high-quality viral genomes",
            proviruses_fna:      "CheckV proviruses.fna: proviral extracts with trimmed host flanking regions"
        }
    }

    parameter_meta {
        virus_fna:             "Viral sequences FASTA from an upstream detection tool (geNomad or VirSorter2)"
        checkv_db_tgz:         "CheckV database as a single compressed archive (.tar.gz or .tar.zst); extracted at runtime"
        sample_name:           "Sample identifier used as output file prefix"
        tool_prefix:           "Short tool identifier prepended to output filenames (e.g. 'genomad' or 'vs2')"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   virus_fna
        File   checkv_db_tgz
        String sample_name
        String tool_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(3.0 * size(checkv_db_tgz, "GB")) + ceil(5.0 * size(virus_fna, "GB"))

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

        # Extract CheckV database
        mkdir -p checkv_db
        if [[ "~{checkv_db_tgz}" == *.tar.zst ]]; then
            zstd -d ~{checkv_db_tgz} --stdout | tar -x -C checkv_db
        else
            tar -xzf ~{checkv_db_tgz} -C checkv_db
        fi
        DB_PATH=$(find checkv_db -name 'genome_db' -type d | head -1 | xargs dirname)
        [[ -z "${DB_PATH}" ]] && DB_PATH=$(ls -d checkv_db/*/ | head -1)

        # Replace spaces with underscores in headers — CheckV truncates at whitespace,
        # causing duplicate IDs when the same base ID appears in multiple tool outputs
        sed 's/ /_/g' ~{virus_fna} > cleaned_input.fna

        # CheckV crashes (Prodigal exits 1) on an empty FASTA — short-circuit cleanly
        N_SEQS=$(grep -c '^>' cleaned_input.fna || true)
        if [[ "${N_SEQS}" -eq 0 ]]; then
            echo "WARNING: empty viral FASTA — skipping CheckV, writing empty outputs" >&2
            touch ~{sample_name}.~{tool_prefix}_checkv_quality.tsv
            touch ~{sample_name}.~{tool_prefix}_checkv_viruses.fna
            touch ~{sample_name}.~{tool_prefix}_checkv_proviruses.fna
            exit 0
        fi

        checkv end_to_end \
            cleaned_input.fna \
            checkv_out \
            -t "${NUM_CPUS}" \
            -d "${DB_PATH}"

        cp checkv_out/quality_summary.tsv ~{sample_name}.~{tool_prefix}_checkv_quality.tsv
        cp checkv_out/viruses.fna         ~{sample_name}.~{tool_prefix}_checkv_viruses.fna
        cp checkv_out/proviruses.fna      ~{sample_name}.~{tool_prefix}_checkv_proviruses.fna
    >>>

    output {
        File quality_summary_tsv = "~{sample_name}.~{tool_prefix}_checkv_quality.tsv"
        File viruses_fna         = "~{sample_name}.~{tool_prefix}_checkv_viruses.fna"
        File proviruses_fna      = "~{sample_name}.~{tool_prefix}_checkv_proviruses.fna"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             16,
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
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
