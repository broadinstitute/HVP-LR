version 1.0

import "../../structs/Structs.wdl"

# ============================================================================
# Pyrodigal-gv — Python/Cython wrapper around Prodigal-gv. Calls open reading
# frames on viral nucleotide sequences. Use it on any nucleotide input that
# upstream tools have NOT already produced protein FASTAs for:
#  - VirSorter2 viral_combined_fa (VS2 emits viral nucleotide but no proteins)
#  - full assembly contigs (incl. contigs neither geNomad nor VS2 flagged —
#    the novel-virus discovery surface)
#  - rescued reads (HvpReadRescue output)
#
# geNomad already runs its own ORF caller internally and surfaces
# virus_proteins.faa, so do NOT run this task on geNomad's virus.fna —
# you'll waste compute and duplicate calls.
#
# pyrodigal-gv ships with the viral training models that standard Prodigal
# lacks (NCLDV: Acanthamoeba polyphaga mimivirus, Chlorella viruses;
# VirSorter2's NCLDV model; Topaz/Agate with genetic code 15; gut phage
# models). CLI defaults to `-p meta` mode — correct for the metagenomic
# inputs above.
# ============================================================================

task PyrodigalGvCallOrfs {

    meta {
        description: "Call ORFs on a viral nucleotide FASTA using pyrodigal-gv (Prodigal-gv via Cython, defaults to metagenomic mode with the bundled viral training models). Emits protein FASTA (for downstream foldseek / mmseqs2 search), nucleotide CDS FASTA, GFF coordinates, and the predicted-gene count."

        tool:          "pyrodigal-gv"
        tool_version:  "0.3.2"
        tool_url:      "https://github.com/althonos/pyrodigal-gv"
        tool_citation: "Larralde M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software. 2022;7(72):4296."

        outputs: {
            proteins_faa: "Translated protein sequences for all called ORFs (FASTA, one entry per gene).",
            cds_fna:      "Nucleotide CDS sequences for all called ORFs (FASTA).",
            genes_gff:    "Per-gene coordinates and metadata (GFF3).",
            num_genes:    "Total number of predicted genes (count of `>` headers in proteins_faa)."
        }
    }

    parameter_meta {
        nucleotide_fasta:      "Input nucleotide FASTA. May be .gz; the task will decompress in place. Typical sources: VirSorter2 viral_combined_fa, full assembly contigs, rescued reads."
        sample_name:           "Sample identifier used as output file prefix"
        source_label:          "Short label identifying the upstream tool / data source (e.g. 'vs2', 'assembly', 'rescued'). Appended to the sample prefix to keep outputs distinguishable when the same sample's nucleotides are called from multiple upstream sources in the same workflow run."
        genetic_code:          "Override genetic code (NCBI translation table). Default 0 = let pyrodigal-gv auto-detect per contig from its bundled training models (recommended for viral input — code 4, 11, 15 etc. vary across viral lineages)."
        closed_ends:           "Pass -c to disallow genes running off contig edges. Default false (allow edge genes) — better for fragmented metagenomic assemblies."
        no_stop_codon:         "Pass --no-stop-codon to suppress '*' at the end of translated proteins. Default false (keep '*' so downstream consumers can distinguish complete genes from partial)."
        extra_args:            "Additional command-line args appended verbatim to the pyrodigal-gv invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   nucleotide_fasta
        String sample_name
        String source_label = "pgv"

        Int     genetic_code  = 0
        Boolean closed_ends   = false
        Boolean no_stop_codon = false

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(15.0 * size(nucleotide_fasta, "GB"))

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

        # Decompress gzipped input — pyrodigal-gv handles both plain and gz
        # transparently in newer releases, but force-decompress here so the
        # output path is uniform and pre-3.x callers behave identically.
        if [[ "~{nucleotide_fasta}" == *.gz ]]; then
            gunzip -c ~{nucleotide_fasta} > input.fa
        else
            cp ~{nucleotide_fasta} input.fa
        fi

        # Build the flag list. Conditional flags emitted only when true so the
        # CLI line stays clean when defaults are used.
        FLAGS=()
        if [[ "~{genetic_code}" -ne 0 ]]; then
            FLAGS+=(-g "~{genetic_code}")
        fi
        if [[ "~{closed_ends}" == "true" ]]; then
            FLAGS+=(-c)
        fi
        if [[ "~{no_stop_codon}" == "true" ]]; then
            FLAGS+=(--no-stop-codon)
        fi

        pyrodigal-gv \
            -i input.fa \
            -a ~{sample_name}.~{source_label}.proteins.faa \
            -d ~{sample_name}.~{source_label}.cds.fna \
            -o ~{sample_name}.~{source_label}.genes.gff \
            -f gff \
            -j "${NUM_CPUS}" \
            "${FLAGS[@]}" \
            ~{extra_args}

        # Gene count — header-count of the protein FASTA is the canonical
        # gene tally. `grep -c` would exit 1 on a zero-gene input and trip
        # set -e; use awk so empty output gives 0.
        awk '/^>/{n++} END{print n+0}' ~{sample_name}.~{source_label}.proteins.faa > num_genes.txt
    >>>

    output {
        File proteins_faa = "~{sample_name}.~{source_label}.proteins.faa"
        File cds_fna      = "~{sample_name}.~{source_label}.cds.fna"
        File genes_gff    = "~{sample_name}.~{source_label}.genes.gff"
        Int  num_genes    = read_int("num_genes.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/pyrodigal-gv:0.3.2"
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
