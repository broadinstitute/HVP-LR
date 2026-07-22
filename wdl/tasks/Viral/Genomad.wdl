version 1.0

import "../../structs/Structs.wdl"

task Genomad {

    meta {
        description: "Run geNomad end-to-end on a merged assembly+rescue FASTA to identify viral sequences and plasmids. Input FASTA is symlinked to a fixed name (input_contigs.fasta.gz) so output subdirectory names are predictable regardless of the upstream filename. Outputs the virus summary TSV (scores, topology, taxonomy) and the viral sequence FASTA for downstream CheckV quality assessment."

        tool:         "geNomad"
        tool_version: "1.9"
        tool_url:     "https://github.com/apcamargo/genomad"
        tool_citation: "Camargo AP, et al. Identification of mobile genetic elements with geNomad. Nat Biotechnol. 2023."

        outputs: {
            virus_summary_tsv:  "geNomad virus summary TSV (one row per viral sequence: score, FDR, hallmarks, topology, taxonomy)",
            virus_fna:          "Viral nucleotide sequences FASTA produced by geNomad; input to CheckV",
            virus_proteins_faa: "Predicted protein sequences (amino-acid FASTA) for ORFs on viral contigs; ready for foldseek / structural search via ProstT5",
            virus_genes_tsv:    "Per-gene TSV from geNomad annotate (gene coordinates, strand, marker hits, annotations)"
        }
    }

    parameter_meta {
        merged_fa_gz:          "Merged assembly+rescue FASTA (gzipped) from HvpReadRescue"
        genomad_db_tgz:        "geNomad database as a single compressed archive (.tar.gz or .tar.zst); extracted at runtime"
        sample_name:           "Sample identifier used as output file prefix"
        extra_args:            "Additional command-line args appended verbatim to the genomad invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   merged_fa_gz
        File   genomad_db_tgz
        String sample_name

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(genomad_db_tgz, "GB")) + ceil(10.0 * size(merged_fa_gz, "GB"))

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

        # Extract geNomad database
        mkdir -p genomad_db
        if [[ "~{genomad_db_tgz}" == *.tar.zst ]]; then
            zstd -d ~{genomad_db_tgz} --stdout | tar -x -C genomad_db
        else
            tar -xzf ~{genomad_db_tgz} -C genomad_db
        fi
        # version.txt is unique to the geNomad DB root — avoids matching the
        # extraction directory itself (also named 'genomad_db').
        DB_PATH=$(find genomad_db -name 'version.txt' -type f | head -1 | xargs -r dirname)
        [[ -z "${DB_PATH}" ]] && DB_PATH=$(ls -d genomad_db/*/ 2>/dev/null | head -1)

        # Symlink input to a fixed name so output subdirectory names are predictable
        ln -sf ~{merged_fa_gz} input_contigs.fasta.gz

        genomad end-to-end \
            --threads "${NUM_CPUS}" \
            ~{extra_args} \
            input_contigs.fasta.gz \
            genomad_out \
            "${DB_PATH}"

        # Rename outputs with sample prefix.
        # The virus FASTA is absent when geNomad finds no viral sequences; touch an
        # empty file so WDL output declaration is always satisfied.
        SUMMARY_TSV=$(find genomad_out -name '*_virus_summary.tsv' 2>/dev/null | head -1)
        if [[ -n "${SUMMARY_TSV}" ]]; then
            cp "${SUMMARY_TSV}" ~{sample_name}.genomad_virus_summary.tsv
        else
            touch ~{sample_name}.genomad_virus_summary.tsv
            echo "WARNING: geNomad produced no virus_summary.tsv — no viral sequences found" >&2
        fi

        VIRUS_FNA=$(find genomad_out -name '*_virus.fna' 2>/dev/null | head -1)
        if [[ -n "${VIRUS_FNA}" ]]; then
            cp "${VIRUS_FNA}" ~{sample_name}.genomad_virus.fna
        else
            touch ~{sample_name}.genomad_virus.fna
            echo "WARNING: geNomad produced no virus.fna — no viral sequences found" >&2
        fi

        # Protein FASTA from geNomad's annotate module. Always produced when
        # the input has any predicted ORFs, even if no viral contigs are
        # ultimately called. Empty-touch fallback keeps the output declaration
        # satisfied for degenerate inputs.
        VIRUS_PROTEINS_FAA=$(find genomad_out -name '*_virus_proteins.faa' 2>/dev/null | head -1)
        if [[ -n "${VIRUS_PROTEINS_FAA}" ]]; then
            cp "${VIRUS_PROTEINS_FAA}" ~{sample_name}.genomad_virus_proteins.faa
        else
            touch ~{sample_name}.genomad_virus_proteins.faa
            echo "WARNING: geNomad produced no virus_proteins.faa — no viral ORFs predicted" >&2
        fi

        VIRUS_GENES_TSV=$(find genomad_out -name '*_virus_genes.tsv' 2>/dev/null | head -1)
        if [[ -n "${VIRUS_GENES_TSV}" ]]; then
            cp "${VIRUS_GENES_TSV}" ~{sample_name}.genomad_virus_genes.tsv
        else
            touch ~{sample_name}.genomad_virus_genes.tsv
            echo "WARNING: geNomad produced no virus_genes.tsv — no viral ORFs predicted" >&2
        fi
    >>>

    output {
        File virus_summary_tsv  = "~{sample_name}.genomad_virus_summary.tsv"
        File virus_fna          = "~{sample_name}.genomad_virus.fna"
        File virus_proteins_faa = "~{sample_name}.genomad_virus_proteins.faa"
        File virus_genes_tsv    = "~{sample_name}.genomad_virus_genes.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
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
