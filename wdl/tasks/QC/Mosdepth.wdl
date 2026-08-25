version 1.0

import "../../structs/Structs.wdl"

task Mosdepth {

    meta {
        description: "Compute windowed sequencing depth from a coordinate-sorted BAM with mosdepth. Emits per-window regions BED (--by), the per-chromosome summary, and the global distribution. Used upstream of the sex-chromosome karyotype classifier (defaults: 1 Mb windows, MAPQ>=20 — the granularity/quality the classifier is calibrated for)."

        tool:          "mosdepth"
        tool_version:  "0.3.11"
        tool_url:      "https://github.com/brentp/mosdepth"
        tool_citation: "Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867-868."

        outputs: {
            regions_bed:     "Per-window mean-depth BED (bgzipped): chrom, start, end, mean depth",
            regions_bed_csi: "CSI index for regions_bed",
            summary_txt:     "Per-chromosome mean-depth summary (mosdepth.summary.txt)",
            global_dist_txt: "Cumulative coverage distribution (mosdepth.global.dist.txt)"
        }
    }

    parameter_meta {
        aligned_bam:           "Coordinate-sorted BAM (aligned to T2T-CHM13v2.0 for the karyotype use)"
        aligned_bai:           "BAI index for aligned_bam"
        prefix:                "Basename for mosdepth outputs"
        window_size:           "Window size in bp for --by (default 1000000)"
        min_mapq:              "Minimum MAPQ to count a read, mosdepth -Q (default 20)"
        extra_args:            "Additional command-line args appended verbatim to the mosdepth invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File aligned_bam
        File aligned_bai
        String prefix

        Int window_size = 1000000
        Int min_mapq = 20

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(aligned_bam, "GB") + size(aligned_bai, "GB")))

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        # Reserve 1 GB for OS + container overhead.
        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        # mosdepth requires the index beside the BAM under a matching name.
        ln -s ~{aligned_bam} reads.bam
        ln -s ~{aligned_bai} reads.bam.bai

        # -n: skip per-base output (only need windowed regions); --by: window size; -Q: MAPQ floor.
        mosdepth -t "${NUM_CPUS}" -n --by ~{window_size} -Q ~{min_mapq} ~{extra_args} ~{prefix} reads.bam
    >>>

    output {
        File regions_bed     = "~{prefix}.regions.bed.gz"
        File regions_bed_csi = "~{prefix}.regions.bed.gz.csi"
        File summary_txt     = "~{prefix}.mosdepth.summary.txt"
        File global_dist_txt = "~{prefix}.mosdepth.global.dist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_1"
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
