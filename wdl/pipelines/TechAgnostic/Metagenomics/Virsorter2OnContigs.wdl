version 1.0

import "../../../tasks/Viral/VirSorter2.wdl" as VS2Tasks

workflow Virsorter2OnContigs {

    meta {
        description: "Run VirSorter2 on a pre-assembled contigs FASTA and emit the viral-combined FASTA needed downstream by HvpViralProteinAnnotation. Tech-agnostic: contigs may come from any short-read or long-read assembler (SPAdes, myloasm, metaFlye, etc.); VS2 operates on sequences, not raw reads. Input FASTA may be plain or gzipped — the underlying VS2 task auto-detects. Only dsDNAphage and ssDNA groups are scanned (matches the long-read HvpViralPipeline configuration)."

        allowNestedInputs: true

        outputs: {
            vs2_score_tsv:         "VirSorter2 final-viral-score.tsv (per-sequence viral scores and group assignments)",
            vs2_viral_combined_fa: "Viral sequences FASTA from VirSorter2 (final-viral-combined.fa); feed directly into HvpViralProteinAnnotation.vs2_viral_combined_fa"
        }
    }

    parameter_meta {
        contigs_fa:      "Assembly contigs FASTA (plain or gzipped). For short-read samples in the HVP cohort this is typically the SPAdes assembly written to the entity's contigs_fasta column."
        sample_name:     "Sample identifier used as output file prefix (e.g. HVP-0007.1_7)"
        vs2_db_tgz:      "VirSorter2 database as a single compressed archive (.tar.gz or .tar.zst)"
        vs2_min_length:  "Minimum sequence length for VirSorter2 (default 1000)"
        vs2_extra_args:  "Additional command-line args appended verbatim to the virsorter run invocation"
    }

    input {
        File   contigs_fa
        String sample_name

        File   vs2_db_tgz

        Int    vs2_min_length = 1000
        String vs2_extra_args = ""
    }

    call VS2Tasks.VirSorter2 as t_01_VirSorter2 {
        input:
            merged_fa_gz = contigs_fa,
            vs2_db_tgz   = vs2_db_tgz,
            sample_name  = sample_name,
            min_length   = vs2_min_length,
            extra_args   = vs2_extra_args
    }

    output {
        File vs2_score_tsv         = t_01_VirSorter2.score_tsv
        File vs2_viral_combined_fa = t_01_VirSorter2.viral_combined_fa
    }
}
