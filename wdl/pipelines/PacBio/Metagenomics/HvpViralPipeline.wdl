version 1.0

import "../../../tasks/Viral/Genomad.wdl"      as GnomadTasks
import "../../../tasks/Viral/VirSorter2.wdl"   as VS2Tasks
import "../../../tasks/Viral/CheckV.wdl"        as CheckVTasks
import "../../../tasks/Viral/ViralSummary.wdl"  as SumTasks

workflow HvpViralPipeline {

    meta {
        description: "Single-sample HVP viral identification and quality assessment pipeline. Accepts the merged assembly+rescue FASTA from HvpReadRescue and runs geNomad and VirSorter2 in parallel to identify viral sequences, then assesses quality of each tool's output with CheckV. Produces a per-contig summary joining all four tool outputs and aggregate scalar metrics for the Terra data table."

        allowNestedInputs: true

        outputs: {
            genomad_virus_summary_tsv:     "geNomad virus_summary.tsv (per-sequence scores, topology, taxonomy)",
            genomad_virus_fna:             "Viral nucleotide sequences FASTA from geNomad; input to CheckV",
            genomad_virus_proteins_faa:    "Predicted protein sequences (amino-acid FASTA) for ORFs on geNomad viral contigs; ready for foldseek structural search via ProstT5",
            genomad_virus_genes_tsv:       "Per-gene TSV from geNomad annotate: gene coordinates, strand, marker hits, annotations",
            vs2_score_tsv:                 "VirSorter2 final-viral-score.tsv (per-sequence viral scores and group assignments)",
            vs2_viral_combined_fa:         "Viral sequences FASTA from VirSorter2; input to CheckV",
            genomad_checkv_quality_tsv:    "CheckV quality_summary.tsv for geNomad viral sequences",
            genomad_checkv_viruses_fna:    "CheckV viruses.fna for geNomad viral sequences (complete/HQ genomes)",
            genomad_checkv_proviruses_fna: "CheckV proviruses.fna for geNomad viral sequences (proviral extracts)",
            vs2_checkv_quality_tsv:        "CheckV quality_summary.tsv for VirSorter2 viral sequences",
            vs2_checkv_viruses_fna:        "CheckV viruses.fna for VirSorter2 viral sequences (complete/HQ genomes)",
            vs2_checkv_proviruses_fna:     "CheckV proviruses.fna for VirSorter2 viral sequences (proviral extracts)",
            viral_contig_summary_tsv:      "Per-contig summary TSV joining geNomad and VirSorter2 scores with CheckV quality results",
            viral_stats_tsv:               "Single-row aggregate viral metrics TSV for the sample",
            n_viral_contigs:               "Total viral contigs called by at least one tool",
            n_hq_viral_contigs:            "Number of high-quality viral contigs (CheckV High-quality or Complete)",
            n_complete_viral:              "Number of complete viral genomes (CheckV Complete)",
            n_both_tools:                  "Number of contigs called viral by both geNomad and VirSorter2",
            n50_viral:                     "N50 length of viral contigs",
            total_viral_bases:             "Total bases across all viral contigs",
            n_1tool:                       "Number of contigs called viral by exactly one tool",
            n_genomad:                     "Number of contigs called viral by geNomad",
            n_vs2:                         "Number of contigs called viral by VirSorter2",
            n_mq:                          "Number of medium-quality viral contigs (CheckV Medium-quality)",
            n_lq:                          "Number of low-quality viral contigs (CheckV Low-quality)",
            n_nd:                          "Number of viral contigs with CheckV quality Not-determined",
            n_no_checkv:                   "Number of viral contigs with no CheckV quality result",
            mean_contig_length:            "Mean contig length across all viral contigs",
            n_dtr:                         "Number of viral contigs with DTR topology",
            n_itr:                         "Number of viral contigs with ITR topology",
            n_provirus:                    "Number of viral contigs classified as provirus",
            n_no_repeats:                  "Number of viral contigs with no terminal repeats",
            mean_completeness:             "Mean CheckV completeness across viral contigs with a completeness estimate",
            mean_contamination:            "Mean CheckV contamination across viral contigs with a contamination estimate"
        }
    }

    parameter_meta {
        merged_fa_gz:       "Merged assembly+rescue FASTA (gzipped) from HvpReadRescue"
        sample_name:        "Sample identifier used as output file prefix (e.g. bc2097)"
        genomad_db_tgz:     "geNomad database as a single compressed archive (.tar.gz or .tar.zst)"
        vs2_db_tgz:         "VirSorter2 database as a single compressed archive (.tar.gz or .tar.zst)"
        checkv_db_tgz:      "CheckV database as a single compressed archive (.tar.gz or .tar.zst)"
        genomad_extra_args: "Additional command-line args appended verbatim to the genomad invocation"
        vs2_min_length:     "Minimum sequence length for VirSorter2 (default 1000)"
        vs2_extra_args:     "Additional command-line args appended verbatim to the virsorter run invocation"
    }

    input {
        File   merged_fa_gz
        String sample_name

        File   genomad_db_tgz
        File   vs2_db_tgz
        File   checkv_db_tgz

        String genomad_extra_args = ""
        Int    vs2_min_length     = 1000
        String vs2_extra_args     = ""
    }

    # geNomad and VirSorter2 run in parallel — both take the merged assembly+rescue FASTA
    call GnomadTasks.Genomad as t_01_Genomad {
        input:
            merged_fa_gz   = merged_fa_gz,
            genomad_db_tgz = genomad_db_tgz,
            sample_name    = sample_name,
            extra_args     = genomad_extra_args
    }

    call VS2Tasks.VirSorter2 as t_02_VirSorter2 {
        input:
            merged_fa_gz = merged_fa_gz,
            vs2_db_tgz   = vs2_db_tgz,
            sample_name  = sample_name,
            min_length   = vs2_min_length,
            extra_args   = vs2_extra_args
    }

    # CheckV runs once per upstream tool output
    call CheckVTasks.CheckV as t_03_CheckVGenomd {
        input:
            virus_fna     = t_01_Genomad.virus_fna,
            checkv_db_tgz = checkv_db_tgz,
            sample_name   = sample_name,
            tool_prefix   = "genomad"
    }

    call CheckVTasks.CheckV as t_04_CheckVVS2 {
        input:
            virus_fna     = t_02_VirSorter2.viral_combined_fa,
            checkv_db_tgz = checkv_db_tgz,
            sample_name   = sample_name,
            tool_prefix   = "vs2"
    }

    # Join all four tool outputs into a per-contig summary
    call SumTasks.ViralContigSummary as t_05_ViralContigSummary {
        input:
            genomad_summary_tsv = t_01_Genomad.virus_summary_tsv,
            genomad_checkv_tsv  = t_03_CheckVGenomd.quality_summary_tsv,
            vs2_score_tsv       = t_02_VirSorter2.score_tsv,
            vs2_checkv_tsv      = t_04_CheckVVS2.quality_summary_tsv,
            sample_name         = sample_name
    }

    # Aggregate per-contig summary to scalar metrics for the data table
    call SumTasks.ViralOverallSummary as t_06_ViralOverallSummary {
        input:
            viral_contig_summary_tsv = t_05_ViralContigSummary.viral_contig_summary_tsv,
            sample_name              = sample_name
    }

    output {
        File  genomad_virus_summary_tsv     = t_01_Genomad.virus_summary_tsv
        File  genomad_virus_fna             = t_01_Genomad.virus_fna
        File  genomad_virus_proteins_faa    = t_01_Genomad.virus_proteins_faa
        File  genomad_virus_genes_tsv       = t_01_Genomad.virus_genes_tsv

        File  vs2_score_tsv                 = t_02_VirSorter2.score_tsv
        File  vs2_viral_combined_fa         = t_02_VirSorter2.viral_combined_fa

        File  genomad_checkv_quality_tsv    = t_03_CheckVGenomd.quality_summary_tsv
        File  genomad_checkv_viruses_fna    = t_03_CheckVGenomd.viruses_fna
        File  genomad_checkv_proviruses_fna = t_03_CheckVGenomd.proviruses_fna

        File  vs2_checkv_quality_tsv        = t_04_CheckVVS2.quality_summary_tsv
        File  vs2_checkv_viruses_fna        = t_04_CheckVVS2.viruses_fna
        File  vs2_checkv_proviruses_fna     = t_04_CheckVVS2.proviruses_fna

        File  viral_contig_summary_tsv      = t_05_ViralContigSummary.viral_contig_summary_tsv

        File  viral_stats_tsv               = t_06_ViralOverallSummary.viral_stats_tsv
        Int   n_viral_contigs               = t_06_ViralOverallSummary.n_viral_contigs
        Int   n_hq_viral_contigs            = t_06_ViralOverallSummary.n_hq_viral
        Int   n_complete_viral              = t_06_ViralOverallSummary.n_complete_viral
        Int   n_both_tools                  = t_06_ViralOverallSummary.n_both_tools
        Int   n50_viral                     = t_06_ViralOverallSummary.n50_viral
        Int   total_viral_bases             = t_06_ViralOverallSummary.total_viral_bases
        Int   n_1tool                       = t_06_ViralOverallSummary.n_1tool
        Int   n_genomad                     = t_06_ViralOverallSummary.n_genomad
        Int   n_vs2                         = t_06_ViralOverallSummary.n_vs2
        Int   n_mq                          = t_06_ViralOverallSummary.n_mq
        Int   n_lq                          = t_06_ViralOverallSummary.n_lq
        Int   n_nd                          = t_06_ViralOverallSummary.n_nd
        Int   n_no_checkv                   = t_06_ViralOverallSummary.n_no_checkv
        Float mean_contig_length            = t_06_ViralOverallSummary.mean_contig_length
        Int   n_dtr                         = t_06_ViralOverallSummary.n_dtr
        Int   n_itr                         = t_06_ViralOverallSummary.n_itr
        Int   n_provirus                    = t_06_ViralOverallSummary.n_provirus
        Int   n_no_repeats                  = t_06_ViralOverallSummary.n_no_repeats
        Float mean_completeness             = t_06_ViralOverallSummary.mean_completeness
        Float mean_contamination            = t_06_ViralOverallSummary.mean_contamination
    }
}
