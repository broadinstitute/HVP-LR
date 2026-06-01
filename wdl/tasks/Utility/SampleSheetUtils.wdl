version 1.0

import "../../structs/Structs.wdl"

task ParseSampleSheet {

    meta {
        description: "Parse a sample sheet TSV (with header) and extract columns as parallel arrays."

        outputs: {
            bam_paths: "Per-sample BAM paths (column 1)",
            genera:    "Per-sample genus names (column 2)",
            species:   "Per-sample species names (column 3)",
            strains:   "Per-sample strain identifiers (column 4)",
            barcodes:  "Per-sample barcode identifiers parsed from BAM basename (bc<digits>)"
        }
    }

    parameter_meta {
        sample_tsv:            "Input TSV with columns: bam, genus, species, strain. For metagenomic samples genus and species may be empty."
        extra_args:            "Additional command-line args appended verbatim to the cut invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File sample_tsv

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(sample_tsv, "GB"))

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

        tail -n +2 ~{sample_tsv} | cut ~{extra_args} -f1 > bams.txt
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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

task MergeSampleStats {

    meta {
        description: "Merge per-sample HifiReadStats TSVs into a single multi-row TSV, adding bam and strain columns from the original sample sheet."

        outputs: {
            merged_stats_tsv: "Concatenated all-samples TSV with header from first input plus appended bam and strain columns"
        }
    }

    parameter_meta {
        sample_tsvs:           "Per-sample HifiReadStats TSV files"
        bam_paths:             "Per-sample BAM paths to append as a column"
        strains:               "Per-sample strain identifiers to append as a column"
        extra_args:            "Additional command-line args appended verbatim to the paste invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[File]   sample_tsvs
        Array[String] bam_paths
        Array[String] strains

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    File bam_paths_file = write_lines(bam_paths)
    File strains_file   = write_lines(strains)

    Int disk_size = 10 + ceil(5.0 * size(sample_tsvs, "GB"))

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

        # Build merged TSV body (header from first file, data rows from all)
        head -n 1 ~{sample_tsvs[0]} > merged_body.tsv
        for f in ~{sep=' ' sample_tsvs}; do
            tail -n +2 "$f" >> merged_body.tsv
        done

        # Build extra columns file (bam, strain)
        printf "bam\tstrain\n" > extra_cols.tsv
        paste ~{extra_args} ~{bam_paths_file} ~{strains_file} >> extra_cols.tsv

        # Combine
        paste merged_body.tsv extra_cols.tsv > all_samples.hifi_read_stats.tsv
    >>>

    output {
        File merged_stats_tsv = "all_samples.hifi_read_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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
