version 1.0

import "../../../tasks/Utility/Taxonomy.wdl" as Tax
import "../../../tasks/Utility/PacBioPbi.wdl" as Pbi
import "../../../tasks/Preprocessing/BamConversion.wdl" as Prep
import "../../../tasks/QC/HifiQC.wdl" as QC

workflow HifiPipeline {

    meta {
        description: "Single-sample PacBio HiFi read QC pipeline. Takes one BAM plus genus/species/sample_name, performs taxonomy lookup, converts BAM to FASTQ, computes seqkit and kraken2 stats, and aggregates per-sample metrics. Returns each metric as an individual scalar plus a human-readable report, the FASTQ, and the raw kraken2 report."

        allowNestedInputs: true

        outputs: {
            sample_fastq:             "Gzipped FASTQ produced from the input BAM",
            sample_kraken_report:     "Raw kraken2 report",
            sample_stats_report:      "Per-sample human-readable formatted report",
            pbi_stats_map:            "Two-column TSV (key<TAB>value) of all PBI-derived metrics (ported from broadinstitute/long-read-pipelines)",
            pbi_reads:                "Number of HiFi reads at or above the PBI quality threshold",
            pbi_bases:                "Total bases across all reads (from PBI qStart/qEnd)",
            pbi_mean_qual:            "Mean per-read Phred quality score from the PBI readQual field",
            pbi_median_qual:          "Median per-read Phred quality score from the PBI readQual field",
            polymerase_mean:          "Mean polymerase read length (sum of subread lengths per ZMW)",
            polymerase_median:        "Median polymerase read length",
            polymerase_stdev:         "Standard deviation of polymerase read lengths",
            polymerase_n50:           "N50 polymerase read length",
            subread_mean:             "Mean subread length (one entry per PBI record)",
            subread_median:           "Median subread length",
            subread_stdev:            "Standard deviation of subread lengths",
            subread_n50:              "N50 subread length",
            tax_id:                   "NCBI taxonomy ID resolved from (genus, species)",
            expected_genome_size:     "Expected genome size in bases from the NCBI species_genome_size table; 'NA' if no entry",
            num_reads:                "Total number of HiFi reads",
            bases_in_reads:           "Total bases across all reads",
            estimate_cvg:             "Estimated coverage = bases_in_reads / expected_genome_size; '-' if expected_genome_size is not available",
            bases_in_reads_over_10kb: "Sum of read lengths for reads longer than 10kb",
            estimate_cvg_reads_10kb:  "Estimated coverage from reads >10kb; '-' if expected_genome_size is not available",
            bases_in_reads_over_20kb: "Sum of read lengths for reads longer than 20kb",
            estimate_cvg_reads_20kb:  "Estimated coverage from reads >20kb; '-' if expected_genome_size is not available",
            q1_read_length:           "Q1 of read length distribution",
            median_read_length:       "Median read length",
            q3_read_length:           "Q3 of read length distribution",
            n50_read_length:          "N50 read length",
            max_read_length:          "Maximum read length",
            pct_q20_bases:            "Percentage of bases at Q20 or higher",
            pct_q30_bases:            "Percentage of bases at Q30 or higher",
            mean_read_accuracy:       "Mean per-read accuracy (percentage) from rq:f tags",
            mean_qual_score:          "Mean per-read Phred quality score derived from mean_read_accuracy",
            mean_passes:              "Mean number of CCS subreads per read from np:i tags",
            mean_read_gc:             "Mean GC content (percentage) across reads",
            pct_bacteria:             "Percentage of reads classified as Bacteria",
            pct_fungi:                "Percentage of reads classified as Fungi",
            pct_virus:                "Percentage of reads classified as Viruses",
            pct_human:                "Percentage of reads classified as Homo sapiens",
            pct_unclassified:         "Percentage of reads left unclassified",
            top_genus:                "Genus with the highest kraken2 classification percentage",
            pct_top_genus:            "Percentage of reads assigned to top_genus",
            top_species:              "Species with the highest kraken2 classification percentage",
            pct_top_species:          "Percentage of reads assigned to top_species"
        }
    }

    parameter_meta {
        input_bam:          "PacBio HiFi reads BAM file for the sample"
        input_pbi:          "PacBio .pbi index file for the input BAM (used for ZMW/polymerase metrics)"
        sample_name:        "Sample identifier used as the output prefix for per-sample artifacts (e.g. bc2019)"
        genus:              "Genus name. Empty string allowed for metagenomic samples (disables taxonomy/coverage estimation)."
        species:            "Species name. Empty string allowed for metagenomic samples."
        kraken2_db_hash:    "Kraken2 database hash.k2d file (pre-extracted)"
        kraken2_db_opts:    "Kraken2 database opts.k2d file (pre-extracted)"
        kraken2_db_taxo:    "Kraken2 database taxo.k2d file (pre-extracted)"
        kraken2_confidence: "Kraken2 confidence threshold"
    }

    input {
        File   input_bam
        File   input_pbi
        String sample_name
        String genus
        String species

        File  kraken2_db_hash    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/hash.k2d"  # !FileCoercion
        File  kraken2_db_opts    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/opts.k2d"  # !FileCoercion
        File  kraken2_db_taxo    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/taxo.k2d"  # !FileCoercion
        Float kraken2_confidence = 0.001
    }

    call Tax.GetTaxIdAndGenomeSize as t_01_GetTaxIdAndGenomeSize {
        input:
            genus   = genus,
            species = species
    }

    call Prep.BamToFastqAndStats as t_02_BamToFastqAndStats {
        input:
            input_bam = input_bam
    }

    call QC.HifiSeqkitStats as t_03_HifiSeqkitStats {
        input:
            input_fastq = t_02_BamToFastqAndStats.fastq_gz
    }

    call QC.HifiKraken2 as t_04_HifiKraken2 {
        input:
            input_fastq     = t_02_BamToFastqAndStats.fastq_gz,
            kraken2_db_hash = kraken2_db_hash,
            kraken2_db_opts = kraken2_db_opts,
            kraken2_db_taxo = kraken2_db_taxo,
            confidence      = kraken2_confidence
    }

    call Pbi.SummarizeHifiPbi as t_05_SummarizeHifiPbi {
        input:
            pbi = input_pbi
    }

    call QC.CreateHifiQCReport as t_06_CreateHifiQCReport {
        input:
            sample_name              = sample_name,
            genus                    = genus,
            species                  = species,
            tax_id                   = t_01_GetTaxIdAndGenomeSize.tax_id,
            expected_genome_size     = t_01_GetTaxIdAndGenomeSize.expected_genome_size,

            num_reads                = t_03_HifiSeqkitStats.num_reads,
            bases_in_reads           = t_03_HifiSeqkitStats.bases_in_reads,
            bases_in_reads_over_10kb = t_03_HifiSeqkitStats.bases_in_reads_over_10kb,
            bases_in_reads_over_20kb = t_03_HifiSeqkitStats.bases_in_reads_over_20kb,
            q1_read_length           = t_03_HifiSeqkitStats.q1_read_length,
            median_read_length       = t_03_HifiSeqkitStats.median_read_length,
            q3_read_length           = t_03_HifiSeqkitStats.q3_read_length,
            n50_read_length          = t_03_HifiSeqkitStats.n50_read_length,
            max_read_length          = t_03_HifiSeqkitStats.max_read_length,
            pct_q20_bases            = t_03_HifiSeqkitStats.pct_q20_bases,
            pct_q30_bases            = t_03_HifiSeqkitStats.pct_q30_bases,
            mean_read_gc             = t_03_HifiSeqkitStats.mean_read_gc,

            mean_read_accuracy       = t_02_BamToFastqAndStats.mean_read_accuracy,
            mean_qual_score          = t_02_BamToFastqAndStats.mean_qual_score,
            mean_passes              = t_02_BamToFastqAndStats.mean_passes,

            pct_bacteria             = t_04_HifiKraken2.pct_bacteria,
            pct_fungi                = t_04_HifiKraken2.pct_fungi,
            pct_virus                = t_04_HifiKraken2.pct_virus,
            pct_human                = t_04_HifiKraken2.pct_human,
            pct_unclassified         = t_04_HifiKraken2.pct_unclassified,
            top_genus                = t_04_HifiKraken2.top_genus,
            pct_top_genus            = t_04_HifiKraken2.pct_top_genus,
            top_species              = t_04_HifiKraken2.top_species,
            pct_top_species          = t_04_HifiKraken2.pct_top_species
    }

    output {
        File   sample_fastq             = t_02_BamToFastqAndStats.fastq_gz
        File   sample_kraken_report     = t_04_HifiKraken2.kraken_report
        File   sample_stats_report      = t_06_CreateHifiQCReport.hifi_stats_report
        File   pbi_stats_map            = t_05_SummarizeHifiPbi.pbi_stats_map

        Int    pbi_reads                = t_05_SummarizeHifiPbi.pbi_reads
        Int    pbi_bases                = t_05_SummarizeHifiPbi.pbi_bases
        Float  pbi_mean_qual            = t_05_SummarizeHifiPbi.pbi_mean_qual
        Float  pbi_median_qual          = t_05_SummarizeHifiPbi.pbi_median_qual
        Int    polymerase_mean          = t_05_SummarizeHifiPbi.polymerase_mean
        Int    polymerase_median        = t_05_SummarizeHifiPbi.polymerase_median
        Int    polymerase_stdev         = t_05_SummarizeHifiPbi.polymerase_stdev
        Int    polymerase_n50           = t_05_SummarizeHifiPbi.polymerase_n50
        Int    subread_mean             = t_05_SummarizeHifiPbi.subread_mean
        Int    subread_median           = t_05_SummarizeHifiPbi.subread_median
        Int    subread_stdev            = t_05_SummarizeHifiPbi.subread_stdev
        Int    subread_n50              = t_05_SummarizeHifiPbi.subread_n50

        Int    tax_id                   = t_01_GetTaxIdAndGenomeSize.tax_id
        String expected_genome_size     = t_01_GetTaxIdAndGenomeSize.expected_genome_size

        Int    num_reads                = t_03_HifiSeqkitStats.num_reads
        Int    bases_in_reads           = t_03_HifiSeqkitStats.bases_in_reads
        String estimate_cvg             = t_06_CreateHifiQCReport.estimate_cvg
        Int    bases_in_reads_over_10kb = t_03_HifiSeqkitStats.bases_in_reads_over_10kb
        String estimate_cvg_reads_10kb  = t_06_CreateHifiQCReport.estimate_cvg_reads_10kb
        Int    bases_in_reads_over_20kb = t_03_HifiSeqkitStats.bases_in_reads_over_20kb
        String estimate_cvg_reads_20kb  = t_06_CreateHifiQCReport.estimate_cvg_reads_20kb
        Int    q1_read_length           = t_03_HifiSeqkitStats.q1_read_length
        Int    median_read_length       = t_03_HifiSeqkitStats.median_read_length
        Int    q3_read_length           = t_03_HifiSeqkitStats.q3_read_length
        Int    n50_read_length          = t_03_HifiSeqkitStats.n50_read_length
        Int    max_read_length          = t_03_HifiSeqkitStats.max_read_length
        Float  pct_q20_bases            = t_03_HifiSeqkitStats.pct_q20_bases
        Float  pct_q30_bases            = t_03_HifiSeqkitStats.pct_q30_bases
        Float  mean_read_accuracy       = t_02_BamToFastqAndStats.mean_read_accuracy
        Float  mean_qual_score          = t_02_BamToFastqAndStats.mean_qual_score
        Int    mean_passes              = t_02_BamToFastqAndStats.mean_passes
        Float  mean_read_gc             = t_03_HifiSeqkitStats.mean_read_gc
        Float  pct_bacteria             = t_04_HifiKraken2.pct_bacteria
        Float  pct_fungi                = t_04_HifiKraken2.pct_fungi
        Float  pct_virus                = t_04_HifiKraken2.pct_virus
        Float  pct_human                = t_04_HifiKraken2.pct_human
        Float  pct_unclassified         = t_04_HifiKraken2.pct_unclassified
        String top_genus                = t_04_HifiKraken2.top_genus
        Float  pct_top_genus            = t_04_HifiKraken2.pct_top_genus
        String top_species              = t_04_HifiKraken2.top_species
        Float  pct_top_species          = t_04_HifiKraken2.pct_top_species
    }
}
