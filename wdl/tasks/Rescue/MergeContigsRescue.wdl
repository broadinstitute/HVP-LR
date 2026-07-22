version 1.0

import "../../structs/Structs.wdl"

task MergeContigsRescue {

    meta {
        description: "Merge assembly contigs with rescued reads into a single FASTA for viral discovery. Assembly contig headers are kept as-is. Viral rescue read headers are prefixed with 'rescue_v_' and unclassified rescue read headers with 'rescue_u_', so downstream tools (geNomad, VirSorter2, CheckV) can distinguish rescue reads from assembled contigs in their output. The merged FASTA is the primary input to HvpViralPipeline."

        tool:         "seqkit + coreutils"
        tool_version: "seqkit 2.12.0"
        tool_url:     "https://github.com/shenwei356/seqkit"

        outputs: {
            merged_fa_gz:           "Gzipped FASTA: assembly contigs + rescue_v_* viral reads + rescue_u_* unclassified reads",
            merged_stats_tsv:       "seqkit stats output for the merged FASTA (tab-separated, all columns)",
            num_assembly_contigs:   "Number of sequences from the assembly in the merged FASTA",
            num_rescue_viral:       "Number of kraken-viral rescue sequences in the merged FASTA",
            num_rescue_unclassified: "Number of kraken-unclassified rescue sequences in the merged FASTA",
            total_sequences:        "Total sequences in the merged FASTA (sum of the three counts above)"
        }
    }

    parameter_meta {
        assembly_primary_fa:      "Primary assembly FASTA from HvpAssembly (plain FASTA, uncompressed)"
        rescue_viral_fa_gz:       "Gzipped FASTA of rescued kraken-viral reads from ReadRescue"
        rescue_unclassified_fa_gz: "Gzipped FASTA of rescued kraken-unclassified reads from ReadRescue"
        sample_name:              "Sample identifier used as output file prefix (e.g. bc2097)"
        runtime_attr_override:    "Override the default runtime attributes"
    }

    input {
        File   assembly_primary_fa
        File   rescue_viral_fa_gz
        File   rescue_unclassified_fa_gz
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(10.0 * size(assembly_primary_fa, "GB"))
                       + ceil(5.0  * size(rescue_viral_fa_gz, "GB"))
                       + ceil(5.0  * size(rescue_unclassified_fa_gz, "GB"))

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

        # Count sequences in each input
        awk '/^>/{c++} END{print c+0}' ~{assembly_primary_fa}         > stat.num_assembly_contigs.txt
        zcat ~{rescue_viral_fa_gz}        | awk '/^>/{c++} END{print c+0}' > stat.num_rescue_viral.txt
        zcat ~{rescue_unclassified_fa_gz} | awk '/^>/{c++} END{print c+0}' > stat.num_rescue_unclassified.txt

        # Merge: assembly contigs (headers unchanged) + rescue reads (headers prefixed)
        {
            cat ~{assembly_primary_fa}
            zcat ~{rescue_viral_fa_gz} \
                | awk '/^>/{print ">rescue_v_" substr($0, 2); next} {print}'
            zcat ~{rescue_unclassified_fa_gz} \
                | awk '/^>/{print ">rescue_u_" substr($0, 2); next} {print}'
        } | gzip > ~{sample_name}.merged.fasta.gz

        # Stats on merged FASTA
        seqkit stats -T -a -j "${NUM_CPUS}" ~{sample_name}.merged.fasta.gz \
            > ~{sample_name}.merged_stats.tsv

        # Total sequence count (sum of three parts — avoids parsing seqkit output)
        echo $(( $(cat stat.num_assembly_contigs.txt) \
               + $(cat stat.num_rescue_viral.txt) \
               + $(cat stat.num_rescue_unclassified.txt) )) > stat.total_sequences.txt
    >>>

    output {
        File merged_fa_gz            = "~{sample_name}.merged.fasta.gz"
        File merged_stats_tsv        = "~{sample_name}.merged_stats.tsv"
        Int  num_assembly_contigs    = read_int("stat.num_assembly_contigs.txt")
        Int  num_rescue_viral        = read_int("stat.num_rescue_viral.txt")
        Int  num_rescue_unclassified = read_int("stat.num_rescue_unclassified.txt")
        Int  total_sequences         = read_int("stat.total_sequences.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
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
