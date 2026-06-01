version 1.0

import "../../structs/Structs.wdl"

task BamToFastqAndStats {

    meta {
        description: "Convert a PacBio HiFi BAM file to fastq.gz and parse per-read tags (np:i, rq:f) to compute mean read accuracy, Phred quality score, and mean number of CCS passes."

        tool:          "samtools"
        tool_version:  "1.23"
        tool_url:      "https://www.htslib.org/"
        tool_citation: "Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008."

        outputs: {
            fastq_gz:           "Compressed FASTQ produced by samtools fastq",
            mean_read_accuracy: "Mean read accuracy from rq:f tags, as a percentage (e.g. 99.876)",
            mean_qual_score:    "Mean Phred quality score derived from mean_read_accuracy",
            mean_passes:        "Mean number of CCS subreads per read from np:i tags, rounded down"
        }
    }

    parameter_meta {
        input_bam:             "PacBio HiFi reads BAM file"
        extra_args:            "Additional command-line args appended verbatim to the samtools fastq invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    String bam_basename = basename(input_bam, ".bam")

    Int disk_size = 10 + ceil(5.0 * size(input_bam, "GB"))

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

        samtools fastq \
            -@ "${NUM_CPUS}" \
            ~{extra_args} \
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

        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f1 > mean_read_accuracy.txt
        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f2 > mean_qual_score.txt
        tail -1 ~{bam_basename}.bam_stats.tsv | cut -f3 > mean_passes.txt
    >>>

    output {
        File   fastq_gz           = "~{bam_basename}.fastq.gz"
        Float  mean_read_accuracy = read_float("mean_read_accuracy.txt")
        Float  mean_qual_score    = read_float("mean_qual_score.txt")
        Int    mean_passes        = read_int("mean_passes.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/samtools:1.23--h96c455f_0"
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
