version 1.0

import "../../structs/Structs.wdl"

task RxSexKaryotype {

    meta {
        description: "Call a sex-chromosome karyotype (46,XX / 46,XY / 47,XXY / 47,XYY / 47,XXX / 45,X / 48,XXYY / OTHER) from mosdepth windows via the x_dosage rx_sex classifier. Computes Rx (chrX non-PAR / autosome) and Ry (chrY euchromatin / autosome) dosage and returns the maximum-posterior karyotype with a confidence, on T2T-CHM13v2.0. Off-multiple / mosaic dosages fall to OTHER by design (see the classifier docstring). Input windows must be 1 Mb (see the Mosdepth task)."

        tool:         "rx_sex.py (x_dosage)"
        tool_version: "0.1.0"
        tool_url:     "https://github.com/broadinstitute/HVP-LR"

        outputs: {
            karyotype_tsv:  "Full classifier output TSV (sample, auto_dp, Rx, Rx_CI, Ry, Ry_CI, call, conf, runner_up)",
            karyotype_call: "Maximum-posterior karyotype label (or OTHER / INSUFFICIENT)",
            confidence:     "Posterior probability of the called karyotype",
            rx:             "Rx dosage: median chrX non-PAR depth / median autosome depth",
            ry:             "Ry dosage: median chrY euchromatin depth / median autosome depth",
            auto_depth:     "Median autosomal window depth",
            runner_up:      "Second-ranked karyotype and its posterior (label:prob)"
        }
    }

    parameter_meta {
        regions_bed:           "mosdepth *.regions.bed.gz at 1 Mb windows on T2T-CHM13v2.0"
        sample_name:           "Sample identifier (used as the output prefix)"
        platform:              "Which baked calibrated config to use: short_read | long_read (default short_read). Ignored if config is supplied."
        config:                "Optional TOML config override; defaults to the baked /opt/x_dosage/configs/<platform>.toml"
        extra_args:            "Additional command-line args appended verbatim to the rx_sex.py invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File regions_bed
        String sample_name

        String platform = "short_read"
        File? config

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(regions_bed, "GB"))

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

        # Pick config: explicit override, else the baked per-platform calibrated config.
        CFG="~{if defined(config) then config else '/opt/x_dosage/configs/' + platform + '.toml'}"
        test -f "${CFG}"

        rx_sex.py --config "${CFG}" ~{extra_args} ~{regions_bed} | tee ~{sample_name}.karyotype.tsv

        # Parse the single data row (line 2) into scalar files for typed WDL outputs.
        # Columns: sample auto_dp Rx Rx_CI Ry Ry_CI call conf runner_up
        tail -n +2 ~{sample_name}.karyotype.tsv | head -1 \
          | awk -F'\t' 'NF>=9 {
                print $2 > "auto.txt"; print $3 > "rx.txt"; print $5 > "ry.txt";
                print $7 > "call.txt"; print $8 > "conf.txt"; print $9 > "runnerup.txt"; next
              } {
                print 0 > "auto.txt"; print 0 > "rx.txt"; print 0 > "ry.txt";
                print "INSUFFICIENT" > "call.txt"; print 0 > "conf.txt"; print "NA" > "runnerup.txt"
              }'
    >>>

    output {
        File karyotype_tsv    = "~{sample_name}.karyotype.tsv"
        String karyotype_call = read_string("call.txt")
        Float confidence      = read_float("conf.txt")
        Float rx              = read_float("rx.txt")
        Float ry              = read_float("ry.txt")
        Float auto_depth      = read_float("auto.txt")
        String runner_up      = read_string("runnerup.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/x_dosage:0.1.0"
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
