version 1.0

import "../../structs/Structs.wdl"

# ── Skani ─────────────────────────────────────────────────────────────────────

task Skani {

    meta {
        description: "Search DAS_Tool-refined bins against the GTDB r226 pre-built sketch database using skani to identify closest bacterial/archaeal reference genomes. The database archive contains index.db, markers.bin, and sketches.db at archive root; the taxonomy TSV is a separate input used only by SkaniAnnotate."
        tool:         "skani"
        tool_version: "0.2.2"
        tool_url:     "https://github.com/bluenote-1577/skani"
        tool_citation: "Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023;20(11):1661-1665."
        outputs: {
            results_tsv: "skani search output TSV: per-bin best hits with ANI, alignment fractions, and reference file paths"
        }
    }

    parameter_meta {
        bins:                  "Refined bin FASTA files from DAS_Tool"
        skani_db_tgz:          "Skani GTDB r226 sketch database as tar.zst archive (index.db, markers.bin, sketches.db at archive root)"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        Array[File] bins
        File        skani_db_tgz
        String      sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(3.0 * size(skani_db_tgz, "GB")) + ceil(2.0 * size(bins, "GB"))
    Int mem_size  = ceil(2.0 * size(skani_db_tgz, "GB")) + ceil(4.0 * size(bins, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        # Extract DB — index.db, markers.bin, sketches.db land at skani_db/ root
        mkdir -p skani_db
        if [[ "~{skani_db_tgz}" == *.tar.zst ]]; then
            zstd -d ~{skani_db_tgz} --stdout | tar -x -C skani_db
        else
            tar -xzf ~{skani_db_tgz} -C skani_db
        fi

        # Stage bins to a local directory so we can glob them
        mkdir -p bins_dir
        for f in ~{sep=' ' bins}; do
            ln -sf "$f" "bins_dir/$(basename "$f")"
        done

        bins=( bins_dir/*.fa )
        echo "Searching ${#bins[@]} bins against GTDB r226 sketch DB" >&2

        skani search \
            -q "${bins[@]}" \
            -d skani_db \
            -o ~{sample_name}.skani_results.tsv \
            -t "${NUM_CPUS}"
    >>>

    output {
        File results_tsv = "~{sample_name}.skani_results.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         24,
        mem_gb:            mem_size,
        disk_gb:           disk_size,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
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

# ── SkaniAnnotate ─────────────────────────────────────────────────────────────

task SkaniAnnotate {

    meta {
        description: "Join skani search results with GTDB r226 taxonomy to add a full d__;p__;...;s__ lineage column. Extracts genome accession from Ref_file path and looks up the taxonomy table by both GB_/RS_-prefixed and bare accession keys."
        tool:         "awk (taxonomy join)"
        outputs: {
            annotated_tsv: "skani results TSV with appended gtdb_taxonomy column"
        }
    }

    parameter_meta {
        skani_results_tsv:     "skani search output TSV from the Skani task"
        taxonomy_tsv:          "GTDB r226 combined taxonomy TSV (gtdb_r226_combined_taxonomy.tsv)"
        sample_name:           "Sample identifier used as output file prefix"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File   skani_results_tsv
        File   taxonomy_tsv
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        awk -F'\t' -v OFS='\t' '
        NR == FNR {
            if (NR == 1) next
            accession = $1
            stripped  = accession
            sub(/^(GB_|RS_)/, "", stripped)
            lookup[accession] = $2
            by_acc[stripped]  = $2
            next
        }

        /^Ref_file/ {
            if (!header_done) {
                print $0, "gtdb_taxonomy"
                header_done = 1
            }
            next
        }

        {
            ref = $1
            n = split(ref, parts, "/")
            acc = parts[n]
            sub(/_genomic\.fna\.gz$/, "", acc)
            sub(/\.gz$/, "", acc)

            if      (acc ~ /^GCA_/) full_key = "GB_" acc
            else if (acc ~ /^GCF_/) full_key = "RS_" acc
            else                     full_key = acc

            if      (full_key in lookup) tax = lookup[full_key]
            else if (acc      in by_acc) tax = by_acc[acc]
            else                         tax = "NA"

            print $0, tax
        }
        ' ~{taxonomy_tsv} ~{skani_results_tsv} > ~{sample_name}.skani_annotated.tsv

        n_hits=$(awk 'NR > 1' ~{sample_name}.skani_annotated.tsv | wc -l)
        n_na=$(awk -F'\t' 'NR > 1 && $NF == "NA"' ~{sample_name}.skani_annotated.tsv | wc -l)
        echo "Annotated ${n_hits} rows; ${n_na} with no taxonomy match" >&2
    >>>

    output {
        File annotated_tsv = "~{sample_name}.skani_annotated.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            8,
        disk_gb:           20,
        boot_disk_gb:      25,
        preemptible_tries: 2,
        max_retries:       1,
        docker:            "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
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
