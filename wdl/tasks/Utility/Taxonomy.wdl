version 1.0

import "../../structs/Structs.wdl"

task GetTaxIdAndGenomeSize {

    meta {
        description: "Look up the NCBI taxonomy ID and expected genome size for a single (genus, species) pair. Downloads the NCBI taxonomy dump and species genome size table, resolves the species via taxonkit name2taxid (scientific name first, then unrestricted), and emits a one-row TSV plus scalar values."

        tool:          "taxonkit"
        tool_version:  "0.20.0"
        tool_url:      "https://github.com/shenwei356/taxonkit"
        tool_citation: "Shen W, Ren H. TaxonKit: a practical and efficient NCBI taxonomy toolkit. Journal of Genetics and Genomics. 2021;48(9):844-850."

        outputs: {
            taxid_and_genome_size_tsv: "One-row TSV with header (tax_id, expected_genome_size, genus, species) and a single data row",
            tax_id:                    "NCBI taxonomy ID resolved from (genus, species)",
            expected_genome_size:      "Expected genome size in bases from the NCBI species_genome_size table; 'NA' if no entry exists for the resolved tax_id"
        }
    }

    parameter_meta {
        genus:                 "Genus name"
        species:               "Species name"
        extra_args:            "Additional command-line args appended verbatim to the taxonkit name2taxid invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        String genus
        String species

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

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

        mkdir -p taxdump
        wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf taxdump.tar.gz -C taxdump

        wget -qO species_genome_size.txt \
            "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt"

        species_name="~{genus} ~{species}"

        tax_id=$(echo "$species_name" \
            | taxonkit name2taxid --sci-name --data-dir taxdump ~{extra_args} \
            | cut -f2)

        if [ -z "$tax_id" ]; then
            tax_id=$(echo "$species_name" \
                | taxonkit name2taxid --data-dir taxdump ~{extra_args} \
                | cut -f2)
        fi

        if [ -z "$tax_id" ]; then
            echo "Error: no taxonomy ID found for '$species_name'" >&2
            exit 1
        fi

        expected_genome_size=$(awk -F'\t' -v id="$tax_id" '$1 == id {print $4}' species_genome_size.txt)
        if [ -z "$expected_genome_size" ]; then
            echo "Warning: no genome size entry found for tax_id $tax_id. Setting to NA." >&2
            expected_genome_size="NA"
        fi

        printf "tax_id\texpected_genome_size\tgenus\tspecies\n" > taxid_and_genome_size.tsv
        printf "%s\t%s\t%s\t%s\n" "$tax_id" "$expected_genome_size" "~{genus}" "~{species}" >> taxid_and_genome_size.tsv

        echo "$tax_id" > tax_id.txt
        echo "$expected_genome_size" > expected_genome_size.txt
    >>>

    output {
        File   taxid_and_genome_size_tsv = "taxid_and_genome_size.tsv"
        Int    tax_id                    = read_int("tax_id.txt")
        String expected_genome_size      = read_string("expected_genome_size.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1"
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
