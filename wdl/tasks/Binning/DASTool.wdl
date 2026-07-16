version 1.0

import "../../structs/Structs.wdl"

task DASTool {

    meta {
        description: "Run DAS_Tool to dereplicate and refine bins from MetaBAT2, MaxBin2, and SemiBin2 into a consensus bin set. Converts each binner's Array[File] FASTA inputs to scaffold-to-bin TSVs internally. Requires docopt R package <= 0.7.0 (pinned install — newer versions dropped the 'short' field from Argument, crashing DAS_Tool 1.1.7 on startup)."
        tool:         "DAS_Tool"
        tool_version: "1.1.7"
        tool_url:     "https://github.com/cmks/DAS_Tool"
        tool_citation: "Sieber CMK, et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nat Microbiol. 2018;3(7):836-843."
        outputs: {
            bins:     "Refined bin FASTA files selected by DAS_Tool",
            eval_tsv: "DAS_Tool per-bin score summary (completeness, contamination, size, N50)"
        }
    }

    parameter_meta {
        contigs_fa:            "Assembly contigs FASTA (uncompressed or .gz)"
        metabat2_bins:         "Bin FASTA files from MetaBAT2 (bin.N.fa)"
        maxbin2_bins:          "Bin FASTA files from MaxBin2 (bin.NNN.fasta)"
        semibin2_bins:         "Bin FASTA files from SemiBin2 (SemiBin_N.fa)"
        sample_name:           "Sample identifier used as output file prefix"
        score_threshold:       "Minimum DAS_Tool score for a bin to be retained (default 0.6)"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File        contigs_fa
        Array[File] metabat2_bins
        Array[File] maxbin2_bins
        Array[File] semibin2_bins
        String      sample_name

        Float score_threshold = 0.6

        RuntimeAttr? runtime_attr_override
    }

    Int binner_bins_gb = ceil(size(metabat2_bins, "GB") + size(maxbin2_bins, "GB") + size(semibin2_bins, "GB"))
    Int disk_size = 20 + ceil(3.0 * size(contigs_fa, "GB")) + (3 * binner_bins_gb)

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        # DAS_Tool 1.1.7 requires docopt <= 0.7.0
        Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/docopt/docopt_0.7.0.tar.gz', repos=NULL, type='source')" 2>&1

        if [[ "~{contigs_fa}" == *.gz ]]; then
            gunzip -c ~{contigs_fa} > contigs.fasta
        else
            cp ~{contigs_fa} contigs.fasta
        fi

        # Stage each binner's bins into named directories so basenames are preserved
        mkdir -p metabat2_bins maxbin2_bins semibin2_bins

        for f in ~{sep=' ' metabat2_bins}; do cp "$f" metabat2_bins/; done
        for f in ~{sep=' ' maxbin2_bins};  do cp "$f" maxbin2_bins/;  done
        for f in ~{sep=' ' semibin2_bins}; do cp "$f" semibin2_bins/; done

        # Convert each bin directory to scaffold-to-bin TSV (contig_id <TAB> bin_name)
        fasta_dir_to_s2b() {
            local dir="$1" ext="$2"
            for f in "${dir}"/*.${ext}; do
                [[ -f "$f" ]] || continue
                bin=$(basename "$f" ".${ext}")
                grep '^>' "$f" | sed 's/^>//' | awk -v b="$bin" '{print $1"\t"b}'
            done
        }

        fasta_dir_to_s2b metabat2_bins fa    > metabat2.tsv
        fasta_dir_to_s2b maxbin2_bins  fasta > maxbin2.tsv
        fasta_dir_to_s2b semibin2_bins fa    > semibin2.tsv

        mkdir -p das_tool_out

        DAS_Tool \
            -i metabat2.tsv,maxbin2.tsv,semibin2.tsv \
            -l metabat2,maxbin2,semibin2 \
            -c contigs.fasta \
            -o das_tool_out/das_tool \
            --write_bins \
            --write_bin_evals \
            --score_threshold ~{score_threshold} \
            --threads "${NUM_CPUS}"

        n_bins=$(ls das_tool_out/das_tool_DASTool_bins/*.fa 2>/dev/null | wc -l)
        echo "DAS_Tool produced ${n_bins} refined bins" >&2
    >>>

    output {
        Array[File] bins     = glob("das_tool_out/das_tool_DASTool_bins/*.fa")
        File        eval_tsv = "das_tool_out/das_tool_DASTool_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         8,
        mem_gb:            16,
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
