version 1.0

import "../../structs/Structs.wdl"

task CheckM2 {

    meta {
        description: "Run CheckM2 to assess completeness and contamination of DAS_Tool-refined bins. Accepts an Array[File] of bin FASTAs and stages them into a local directory for CheckM2, which requires a directory input."
        tool:         "CheckM2"
        tool_version: "1.0.2"
        tool_url:     "https://github.com/chklovski/CheckM2"
        tool_citation: "Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods. 2023;20(8):1203-1212."
        outputs: {
            quality_report_tsv: "CheckM2 quality_report.tsv: per-bin completeness, contamination, coding density, contig N50, genome size"
        }
    }

    parameter_meta {
        bins:                  "Refined bin FASTA files from DAS_Tool (extension .fa)"
        checkm2_db_dmnd:       "CheckM2 diamond database (uniref100.KO.1.dmnd) — passed as a direct File, no archive needed"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        Array[File] bins
        File        checkm2_db_dmnd
        String      sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(3.0 * size(checkm2_db_dmnd, "GB")) + ceil(3.0 * size(bins, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        mkdir -p bins_dir
        for f in ~{sep=' ' bins}; do
            ln -sf "$f" "bins_dir/$(basename "$f")"
        done

        checkm2 predict \
            --input bins_dir \
            --output-directory checkm2_out \
            --database_path ~{checkm2_db_dmnd} \
            --extension fa \
            --threads "${NUM_CPUS}"

        cp checkm2_out/quality_report.tsv ~{sample_name}.checkm2_quality.tsv
    >>>

    output {
        File quality_report_tsv = "~{sample_name}.checkm2_quality.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         16,
        mem_gb:            32,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/checkm2:1.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:          select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:       select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb,  default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:   select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:      select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:       select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:           select_first([runtime_attr.docker,            default_attr.docker])
    }
}
