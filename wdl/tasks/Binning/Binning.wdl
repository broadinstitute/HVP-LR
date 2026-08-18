version 1.0

import "../../structs/Structs.wdl"

# ── JgiDepth ──────────────────────────────────────────────────────────────────

task JgiDepth {

    meta {
        description: "Generate per-contig depth/abundance from a sorted BAM using jgi_summarize_bam_contig_depths (MetaBAT2 utility). Output depth TSV is shared by MetaBAT2 and MaxBin2 so the BAM scan happens only once."
        tool:         "MetaBAT2 (jgi_summarize_bam_contig_depths)"
        tool_version: "2.15"
        tool_url:     "https://bitbucket.org/berkeleylab/metabat"
        outputs: {
            depth_txt: "Per-contig depth TSV for use by MetaBAT2 and MaxBin2"
        }
    }

    parameter_meta {
        contigs_fa:            "Assembly contigs FASTA (uncompressed or .gz)"
        sorted_bam:            "Sorted BAM of reads aligned to contigs (from HvpAssembly)"
        sorted_bam_bai:        "BAI index for sorted_bam"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   contigs_fa
        File   sorted_bam
        File   sorted_bam_bai
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(3.0 * size(sorted_bam, "GB")) + ceil(2.0 * size(contigs_fa, "GB"))

    command <<<
        set -euxo pipefail

        # BAM and BAI must be in the same directory
        mkdir -p workdir
        cp ~{sorted_bam}     workdir/align.bam
        cp ~{sorted_bam_bai} workdir/align.bam.bai

        jgi_summarize_bam_contig_depths \
            --outputDepth ~{sample_name}.depth.txt \
            workdir/align.bam
    >>>

    output {
        File depth_txt = "~{sample_name}.depth.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         8,
        mem_gb:            16,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-binning:0.1.0"
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

# ── MetaBAT2 ──────────────────────────────────────────────────────────────────

task MetaBAT2 {

    meta {
        description: "Run MetaBAT2 binning on assembly contigs using pre-computed contig depth from JgiDepth."
        tool:         "MetaBAT2"
        tool_version: "2.15"
        tool_url:     "https://bitbucket.org/berkeleylab/metabat"
        outputs: {
            bins: "Array of bin FASTA files (bin.N.fa); empty array if no bins produced"
        }
    }

    parameter_meta {
        contigs_fa:            "Assembly contigs FASTA (uncompressed or .gz)"
        depth_txt:             "Per-contig depth TSV from JgiDepth"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   contigs_fa
        File   depth_txt
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(contigs_fa, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        if [[ "~{contigs_fa}" == *.gz ]]; then
            gunzip -c ~{contigs_fa} > contigs.fasta
        else
            cp ~{contigs_fa} contigs.fasta
        fi

        mkdir -p metabat2_out

        metabat2 \
            --seed 42 \
            -i contigs.fasta \
            -o metabat2_out/bin \
            -a ~{depth_txt} \
            -t "${NUM_CPUS}"

        n_bins=$(find metabat2_out -maxdepth 1 -name 'bin.*.fa' 2>/dev/null | wc -l)
        echo "MetaBAT2 [~{sample_name}] produced ${n_bins} bins" >&2
    >>>

    output {
        Array[File] bins = glob("metabat2_out/bin.*.fa")
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         8,
        mem_gb:            16,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-binning:0.1.0"
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

# ── MaxBin2 ───────────────────────────────────────────────────────────────────

task MaxBin2 {

    meta {
        description: "Run MaxBin2 binning on assembly contigs using pre-computed contig depth from JgiDepth."
        tool:         "MaxBin2"
        tool_version: "2.2.7"
        tool_url:     "https://sourceforge.net/projects/maxbin2/"
        outputs: {
            bins: "Array of bin FASTA files (bin.NNN.fasta); empty array if no bins produced"
        }
    }

    parameter_meta {
        contigs_fa:            "Assembly contigs FASTA (uncompressed or .gz)"
        depth_txt:             "Per-contig depth TSV from JgiDepth"
        sample_name:           "Sample identifier used as output file prefix"
        min_contig_length:     "Minimum contig length for binning (default 2500)"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   contigs_fa
        File   depth_txt
        String sample_name

        Int min_contig_length = 2500

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(contigs_fa, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        if [[ "~{contigs_fa}" == *.gz ]]; then
            gunzip -c ~{contigs_fa} > contigs.fasta
        else
            cp ~{contigs_fa} contigs.fasta
        fi

        mkdir -p maxbin2_out

        run_MaxBin.pl \
            -contig contigs.fasta \
            -abund ~{depth_txt} \
            -out maxbin2_out/bin \
            -min_contig_length ~{min_contig_length} \
            -thread "${NUM_CPUS}"

        n_bins=$(find maxbin2_out -maxdepth 1 -name 'bin.*.fasta' 2>/dev/null | wc -l)
        echo "MaxBin2 [~{sample_name}] produced ${n_bins} bins" >&2
    >>>

    output {
        Array[File] bins = glob("maxbin2_out/bin.*.fasta")
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         16,
        mem_gb:            32,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-binning:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:          select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:       select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb,  default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:   select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:      select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:       select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:           select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# ── SemiBin2 ──────────────────────────────────────────────────────────────────

task SemiBin2 {

    meta {
        description: "Run SemiBin2 single-sample binning with a pretrained environment model. The semibin_environment string must match one of SemiBin2's built-in environments (human_gut, human_oral, dog_gut, cat_gut, mouse_gut, pig_gut, ocean, soil, built_environment, wastewater, chicken_caecum, global)."
        tool:         "SemiBin2"
        tool_version: "2.1.0"
        tool_url:     "https://github.com/BigDataBiology/SemiBin"
        outputs: {
            bins: "Array of bin FASTA files from SemiBin2 (SemiBin_N.fa)"
        }
    }

    parameter_meta {
        contigs_fa:            "Assembly contigs FASTA (uncompressed or .gz)"
        sorted_bam:            "Sorted BAM of reads aligned to contigs"
        semibin_environment:   "Pretrained SemiBin2 environment model (e.g. human_gut)"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   contigs_fa
        File   sorted_bam
        String semibin_environment
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(contigs_fa, "GB")) + ceil(3.0 * size(sorted_bam, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        if [[ "~{contigs_fa}" == *.gz ]]; then
            gunzip -c ~{contigs_fa} > contigs.fasta
            CONTIGS=contigs.fasta
        else
            CONTIGS=~{contigs_fa}
        fi

        SemiBin2 single_easy_bin \
            -i "${CONTIGS}" \
            -b ~{sorted_bam} \
            -o semibin2_out \
            --environment ~{semibin_environment} \
            --threads "${NUM_CPUS}"

        # SemiBin2 2.2.x writes compressed bins (.fa.gz) — decompress for downstream tools
        for f in semibin2_out/output_bins/SemiBin_*.fa.gz; do
            [[ -f "$f" ]] && gunzip "$f"
        done

        n_bins=$(find semibin2_out/output_bins -maxdepth 1 -name 'SemiBin_*.fa' 2>/dev/null | wc -l)
        echo "SemiBin2 [~{sample_name}] produced ${n_bins} bins" >&2
    >>>

    output {
        Array[File] bins = glob("semibin2_out/output_bins/SemiBin_*.fa")
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         16,
        mem_gb:            32,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-binning:0.1.0"
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
