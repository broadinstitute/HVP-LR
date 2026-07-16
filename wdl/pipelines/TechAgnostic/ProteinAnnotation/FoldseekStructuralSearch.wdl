version 1.0

import "../../../tasks/ProteinAnnotation/Foldseek.wdl" as FS

workflow FoldseekStructuralSearch {

    meta {
        description: "Run foldseek easy-search of a query structure set against a target structure set. Single-task pipeline: stages PDB/mmCIF files into query/ and target/ directories, runs foldseek easy-search end-to-end (createdb, structural prefilter, 3Di + AA alignment, format), and returns the hits TSV plus a scalar hit count. Tech-agnostic; works on any protein structure set regardless of upstream sequencing platform."

        allowNestedInputs: true

        outputs: {
            results_tsv: "Foldseek easy-search alignment hits, one row per (query, target) hit. Columns determined by `format_output`.",
            num_hits:    "Number of rows in results_tsv (count of query-target hit pairs surviving the e-value cutoff)."
        }
    }

    parameter_meta {
        query_structures:  "Per-structure files (PDB or mmCIF, optionally .gz) used as queries"
        target_structures: "Per-structure files (PDB or mmCIF, optionally .gz) used as the search target set"
        prefix:            "Basename for the results TSV (e.g. 'my_run' yields my_run.m8)"
        evalue_cutoff:     "Maximum e-value of reported hits (foldseek -e). Default 0.001."
        format_output:     "Comma-separated foldseek --format-output column spec. Default reports query/target/evalue/bits/fident/alnlen/prob."
    }

    input {
        Array[File] query_structures
        Array[File] target_structures
        String      prefix

        Float  evalue_cutoff = 0.001
        String format_output = "query,target,evalue,bits,fident,alnlen,prob"
    }

    call FS.FoldseekEasySearch as t_01_FoldseekEasySearch {
        input:
            query_structures  = query_structures,
            target_structures = target_structures,
            prefix            = prefix,
            evalue_cutoff     = evalue_cutoff,
            format_output     = format_output
    }

    output {
        File results_tsv = t_01_FoldseekEasySearch.results_tsv
        Int  num_hits    = t_01_FoldseekEasySearch.num_hits
    }
}
