version 1.0

import "../../../tasks/Rescue/ReadRescue.wdl"         as RescueTasks
import "../../../tasks/Rescue/MergeContigsRescue.wdl" as MergeTasks

workflow HvpReadRescue {

    meta {
        description: "Single-sample HVP read rescue workflow. Cross-references per-read kraken2 classifications against the minimap2 read-to-contig BAM to identify reads unaccounted for by the assembly (unmapped or high clip fraction). Partitions those reads into kraken2-viral and kraken2-unclassified sets, then merges them with the assembly contigs into a single FASTA for viral discovery. The merged FASTA (assembly contigs + rescue_v_* viral + rescue_u_* unclassified) is the primary input to HvpViralPipeline."

        allowNestedInputs: true

        outputs: {
            inventory_tsv:              "Per-read TSV: read_id, kraken_status, kraken_taxid, kraken_taxon_name, mapped, mapq, clip_frac, read_len, accounted, category, contig",
            rescue_viral_fa_gz:         "Gzipped FASTA of rescued kraken-viral reads",
            rescue_unclassified_fa_gz:  "Gzipped FASTA of rescued kraken-unclassified reads",
            rescue_summary_tsv:         "Summary TSV: per-category read counts and thresholds",
            merged_fa_gz:               "Merged FASTA: assembly contigs + rescue_v_* + rescue_u_*; primary input to HvpViralPipeline",
            merged_stats_tsv:           "seqkit stats for the merged FASTA",
            num_rescue_viral_merged:    "Number of kraken-viral rescue sequences in the merged FASTA",
            num_rescue_unclassified_merged: "Number of kraken-unclassified rescue sequences in the merged FASTA",
            num_contigs_plus_rescue_reads: "Total sequences in the merged FASTA (assembly contigs + rescue reads)"
        }
    }

    parameter_meta {
        cleaned_bam:           "Cleaned unmapped BAM from HvpReadProcessing (spike-in and human reads removed)"
        sorted_bam:            "Read-to-contig alignment BAM from HvpAssembly (minimap2 map-hifi, coordinate-sorted)"
        sorted_bam_bai:        "BAI index for sorted_bam"
        assembly_primary_fa:   "Primary assembly FASTA from HvpAssembly (myloasm output)"
        kraken_output:         "Per-read kraken2 classification file from HvpReadProcessing"
        kraken_report:         "Kraken2 report from HvpReadProcessing"
        sample_name:           "Sample identifier used as output file prefix (e.g. bc2097)"
        clip_frac_max:         "Maximum (soft+hard) clip fraction for a read to be considered accounted-for by the assembly (default 0.5)"
        min_rescue_read_len:   "Minimum read length in bases for a rescued read to be included in output FASTA (default 1000)"
    }

    input {
        File   cleaned_bam
        File   sorted_bam
        File   sorted_bam_bai
        File   assembly_primary_fa
        File   kraken_output
        File   kraken_report
        String sample_name

        Float  clip_frac_max       = 0.5
        Int    min_rescue_read_len = 1000
    }

    call RescueTasks.ReadRescue as t_01_ReadRescue {
        input:
            sorted_bam           = sorted_bam,
            sorted_bam_bai       = sorted_bam_bai,
            cleaned_bam          = cleaned_bam,
            kraken_output        = kraken_output,
            kraken_report        = kraken_report,
            sample_name          = sample_name,
            clip_frac_max        = clip_frac_max,
            min_rescue_read_len  = min_rescue_read_len
    }

    call MergeTasks.MergeContigsRescue as t_02_MergeContigsRescue {
        input:
            assembly_primary_fa      = assembly_primary_fa,
            rescue_viral_fa_gz       = t_01_ReadRescue.rescue_viral_fa_gz,
            rescue_unclassified_fa_gz = t_01_ReadRescue.rescue_unclassified_fa_gz,
            sample_name              = sample_name
    }

    output {
        File inventory_tsv             = t_01_ReadRescue.inventory_tsv
        File rescue_viral_fa_gz        = t_01_ReadRescue.rescue_viral_fa_gz
        File rescue_unclassified_fa_gz = t_01_ReadRescue.rescue_unclassified_fa_gz
        File rescue_summary_tsv        = t_01_ReadRescue.rescue_summary_tsv

        File merged_fa_gz                      = t_02_MergeContigsRescue.merged_fa_gz
        File merged_stats_tsv                  = t_02_MergeContigsRescue.merged_stats_tsv
        Int  num_rescue_viral_merged           = t_02_MergeContigsRescue.num_rescue_viral
        Int  num_rescue_unclassified_merged    = t_02_MergeContigsRescue.num_rescue_unclassified
        Int  num_contigs_plus_rescue_reads     = t_02_MergeContigsRescue.total_sequences
    }
}
