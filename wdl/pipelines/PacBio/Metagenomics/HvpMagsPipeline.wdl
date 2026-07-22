version 1.0

import "../../../tasks/Binning/Binning.wdl" as BinningTasks
import "../../../tasks/Binning/DASTool.wdl"  as DASToolTasks
import "../../../tasks/MAG/CheckM2.wdl"      as CheckM2Tasks
import "../../../tasks/MAG/Skani.wdl"        as SkaniTasks
import "../../../tasks/MAG/MagSummary.wdl"   as SumTasks

workflow HvpMagsPipeline {

    meta {
        description: "Single-sample HVP metagenome-assembled genome (MAG) pipeline. Accepts the primary assembly FASTA and a sorted BAM from HvpAssembly, bins contigs with MetaBAT2, MaxBin2, and SemiBin2 in parallel, refines bins with DAS_Tool, assesses quality with CheckM2 and skani taxonomy with GTDB r226, and produces per-bin and per-sample aggregate outputs for the Terra data table."

        allowNestedInputs: true

        outputs: {
            dastool_eval_tsv:          "DAS_Tool per-bin score summary (completeness, contamination, size, N50)",
            checkm2_quality_tsv:       "CheckM2 quality_report.tsv: per-bin completeness, contamination, coding density, contig N50, genome size",
            skani_results_tsv:         "skani search output TSV: per-bin best hits with ANI and alignment fractions",
            skani_annotated_tsv:       "skani results TSV with appended GTDB r226 gtdb_taxonomy column",
            bin_summary_tsv:           "Per-bin combined summary: MIMAG quality tier, GTDB taxonomy, CheckM2 and skani metrics",
            mag_stats_tsv:             "Single-row TSV of all assembly and bin aggregate metrics for the sample",
            total_bins:                "Total number of DAS_Tool-refined bins",
            n_hq_mags:                 "Bins with completeness >= 90 and contamination < 5 (MIMAG HQ)",
            n_medium_quality:          "Bins with completeness >= 50 and contamination < 10 (MIMAG MQ)",
            n_low_quality:             "Bins with completeness < 50 and contamination < 10 (MIMAG LQ)",
            n_contaminated:            "Bins with contamination >= 10 (MIMAG contaminated)",
            n_no_skani_hit:            "Bins with no skani reference hit",
            n_single_contig:           "Bins composed of a single contig (potential complete circular MAG)",
            n_single_contig_hq:        "HQ bins composed of a single contig",
            total_bases_all_bins:      "Total base pairs across all bins",
            total_bases_hq:            "Total base pairs in HQ bins",
            total_bases_medium:        "Total base pairs in medium-quality bins",
            pct_assembly_bases_binned: "Percentage of assembly bases captured in any bin",
            pct_assembly_bases_hq:     "Percentage of assembly bases captured in HQ bins",
            mean_completeness_all:     "Mean CheckM2 completeness across all bins",
            mean_completeness_hq:      "Mean CheckM2 completeness across HQ bins",
            mean_contamination_all:    "Mean CheckM2 contamination across all bins",
            mean_contamination_hq:     "Mean CheckM2 contamination across HQ bins",
            mean_quality_score_hq:     "Mean CheckM2 quality score (completeness - 5*contamination) across HQ bins",
            n_distinct_species_all:    "Number of distinct GTDB species-level assignments across all bins",
            n_distinct_species_hq:     "Number of distinct GTDB species-level assignments among HQ bins",
            n_distinct_genera_all:     "Number of distinct GTDB genus-level assignments across all bins",
            n_distinct_genera_hq:      "Number of distinct GTDB genus-level assignments among HQ bins"
        }
    }

    parameter_meta {
        assembly_primary_fa:   "Primary assembly contigs FASTA from HvpAssembly (uncompressed or .gz)"
        sorted_bam:            "Sorted BAM of reads aligned to primary assembly from HvpAssembly"
        sorted_bam_bai:        "BAI index for sorted_bam"
        asm_stats_tsv:         "Assembly seqkit stats TSV from HvpAssembly (this.asm_stats_tsv)"
        sample_name:           "Sample identifier used as output file prefix (e.g. bc2097)"
        semibin_environment:   "Pretrained SemiBin2 environment model (human_gut, human_oral, dog_gut, cat_gut, mouse_gut, pig_gut, ocean, soil, built_environment, wastewater, chicken_caecum, global)"
        skani_db_tgz:          "Skani GTDB r226 sketch database as tar.zst archive"
        skani_taxonomy_tsv:    "GTDB r226 combined taxonomy TSV (gtdb_r226_combined_taxonomy.tsv)"
        checkm2_db_dmnd:       "CheckM2 diamond database (uniref100.KO.1.dmnd)"
        das_tool_score_threshold: "Minimum DAS_Tool score for a bin to be retained (default 0.6)"
    }

    input {
        File   assembly_primary_fa
        File   sorted_bam
        File   sorted_bam_bai
        File   asm_stats_tsv
        String sample_name
        String semibin_environment

        File   skani_db_tgz
        File   skani_taxonomy_tsv
        File   checkm2_db_dmnd

        Float  das_tool_score_threshold = 0.6
    }

    # Compute per-contig depth once — shared by MetaBAT2 and MaxBin2
    call BinningTasks.JgiDepth as t_01_JgiDepth {
        input:
            contigs_fa      = assembly_primary_fa,
            sorted_bam      = sorted_bam,
            sorted_bam_bai  = sorted_bam_bai,
            sample_name     = sample_name
    }

    # Three binners run in parallel using the same depth and contigs
    call BinningTasks.MetaBAT2 as t_02_MetaBAT2 {
        input:
            contigs_fa  = assembly_primary_fa,
            depth_txt   = t_01_JgiDepth.depth_txt,
            sample_name = sample_name
    }

    call BinningTasks.MaxBin2 as t_03_MaxBin2 {
        input:
            contigs_fa  = assembly_primary_fa,
            depth_txt   = t_01_JgiDepth.depth_txt,
            sample_name = sample_name
    }

    call BinningTasks.SemiBin2 as t_04_SemiBin2 {
        input:
            contigs_fa           = assembly_primary_fa,
            sorted_bam           = sorted_bam,
            semibin_environment  = semibin_environment,
            sample_name          = sample_name
    }

    # DAS_Tool dereplicates and refines bins from all three binners
    call DASToolTasks.DASTool as t_05_DASTool {
        input:
            contigs_fa       = assembly_primary_fa,
            metabat2_bins    = t_02_MetaBAT2.bins,
            maxbin2_bins     = t_03_MaxBin2.bins,
            semibin2_bins    = t_04_SemiBin2.bins,
            sample_name      = sample_name,
            score_threshold  = das_tool_score_threshold
    }

    # CheckM2 and skani run in parallel on the DAS_Tool-refined bins
    call CheckM2Tasks.CheckM2 as t_06_CheckM2 {
        input:
            bins            = t_05_DASTool.bins,
            checkm2_db_dmnd = checkm2_db_dmnd,
            sample_name     = sample_name
    }

    call SkaniTasks.Skani as t_07_Skani {
        input:
            bins        = t_05_DASTool.bins,
            skani_db_tgz = skani_db_tgz,
            sample_name = sample_name
    }

    # Annotate skani hits with GTDB r226 taxonomy
    call SkaniTasks.SkaniAnnotate as t_08_SkaniAnnotate {
        input:
            skani_results_tsv = t_07_Skani.results_tsv,
            taxonomy_tsv      = skani_taxonomy_tsv,
            sample_name       = sample_name
    }

    # Per-bin combined summary (CheckM2 + skani taxonomy)
    call SumTasks.BinSummary as t_09_BinSummary {
        input:
            checkm2_quality_tsv = t_06_CheckM2.quality_report_tsv,
            skani_annotated_tsv = t_08_SkaniAnnotate.annotated_tsv,
            sample_name         = sample_name
    }

    # Aggregate to per-sample scalars for the Terra data table
    call SumTasks.MagSummary as t_10_MagSummary {
        input:
            bin_summary_tsv = t_09_BinSummary.bin_summary_tsv,
            asm_stats_tsv   = asm_stats_tsv,
            sample_name     = sample_name
    }

    output {
        File  dastool_eval_tsv          = t_05_DASTool.eval_tsv
        File  checkm2_quality_tsv       = t_06_CheckM2.quality_report_tsv
        File  skani_results_tsv         = t_07_Skani.results_tsv
        File  skani_annotated_tsv       = t_08_SkaniAnnotate.annotated_tsv
        File  bin_summary_tsv           = t_09_BinSummary.bin_summary_tsv
        File  mag_stats_tsv             = t_10_MagSummary.mag_stats_tsv
        Int   total_bins                = t_10_MagSummary.total_bins
        Int   n_hq_mags                 = t_10_MagSummary.n_high_quality
        Int   n_medium_quality          = t_10_MagSummary.n_medium_quality
        Int   n_low_quality             = t_10_MagSummary.n_low_quality
        Int   n_contaminated            = t_10_MagSummary.n_contaminated
        Int   n_no_skani_hit            = t_10_MagSummary.n_no_skani_hit
        Int   n_single_contig           = t_10_MagSummary.n_single_contig
        Int   n_single_contig_hq        = t_10_MagSummary.n_single_contig_hq
        Int   total_bases_all_bins      = t_10_MagSummary.total_bases_all_bins
        Int   total_bases_hq            = t_10_MagSummary.total_bases_hq
        Int   total_bases_medium        = t_10_MagSummary.total_bases_medium
        Float pct_assembly_bases_binned = t_10_MagSummary.pct_assembly_bases_binned
        Float pct_assembly_bases_hq     = t_10_MagSummary.pct_assembly_bases_hq
        Float mean_completeness_all     = t_10_MagSummary.mean_completeness_all
        Float mean_completeness_hq      = t_10_MagSummary.mean_completeness_hq
        Float mean_contamination_all    = t_10_MagSummary.mean_contamination_all
        Float mean_contamination_hq     = t_10_MagSummary.mean_contamination_hq
        Float mean_quality_score_hq     = t_10_MagSummary.mean_quality_score_hq
        Int   n_distinct_species_all    = t_10_MagSummary.n_distinct_species_all
        Int   n_distinct_species_hq     = t_10_MagSummary.n_distinct_species_hq
        Int   n_distinct_genera_all     = t_10_MagSummary.n_distinct_genera_all
        Int   n_distinct_genera_hq      = t_10_MagSummary.n_distinct_genera_hq
    }
}
