version 1.0

import "../../../tasks/Preprocessing/BamConversion.wdl" as Prep
import "../../../tasks/QC/HifiQC.wdl"                   as QC
import "../../../tasks/Preprocessing/ReadCleaning.wdl"  as Clean

workflow HvpReadProcessing {

    meta {
        description: "Single-sample HVP read processing workflow. Converts a raw unmapped PacBio HiFi BAM to FASTQ, computes read QC metrics, removes spike-in reads (minimap2, NEB HiFi spike-in panel), classifies reads with kraken2 (PlusPF DB), and produces a cleaned BAM with spike-in and human reads removed. The cleaned BAM preserves all original PacBio tags (base modifications, rq:f, np:i) and is the primary input to HvpAssembly. FASTQ is used as an internal intermediate format only; the final output currency is BAM."

        allowNestedInputs: true

        outputs: {
            cleaned_bam:                  "Cleaned BAM with spike-in and human reads removed; all PacBio tags preserved; input to HvpAssembly",
            cleaned_bam_bai:              "BAI index for cleaned_bam",
            num_spikein_excluded:         "Number of reads excluded as spike-in",
            num_human_excluded:           "Number of reads excluded as human (taxid 9606)",
            num_reads_cleaned_bam:        "Total reads remaining in cleaned_bam",
            mean_read_accuracy:           "Mean read accuracy from rq:f BAM tags (percentage)",
            mean_qual_score:              "Mean Phred quality score derived from mean_read_accuracy",
            mean_passes:                  "Mean number of CCS subreads per read from np:i BAM tags",
            seqkit_stats_tsv:             "seqkit stats output on raw FASTQ (length stats, N50, Q20/Q30, GC%)",
            seqkit_fx2tab_tsv:            "seqkit fx2tab output: bases in reads >10kb and >20kb",
            num_reads:                    "Total number of HiFi reads",
            bases_in_reads:               "Total bases across all reads",
            max_read_length:              "Maximum read length",
            q1_read_length:               "Q1 (25th percentile) read length",
            median_read_length:           "Median read length",
            q3_read_length:               "Q3 (75th percentile) read length",
            n50_read_length:              "N50 read length",
            pct_q20_bases:                "Percentage of bases at Q20 or higher",
            pct_q30_bases:                "Percentage of bases at Q30 or higher",
            mean_read_gc:                 "Mean GC content (percentage) across reads",
            bases_in_reads_over_10kb:     "Sum of read lengths for reads longer than 10kb",
            bases_in_reads_over_20kb:     "Sum of read lengths for reads longer than 20kb",
            spikein_report:               "Long-format spike-in report: one row per reference with read count and fraction",
            spikein_stats_tsv:            "Wide-format spike-in stats: one column pair per reference",
            total_spikein_reads:          "Total reads identified as spike-in",
            fraction_spikein:             "Fraction of input reads identified as spike-in",
            num_Lambda_NEB_2026_125_bp:       "Reads classified as Lambda phage 125bp spike-in",
            fraction_Lambda_NEB_2026_125_bp:  "Fraction of reads classified as Lambda phage 125bp spike-in",
            num_Lambda_NEB_2026_564_bp:       "Reads classified as Lambda phage 564bp spike-in",
            fraction_Lambda_NEB_2026_564_bp:  "Fraction of reads classified as Lambda phage 564bp spike-in",
            num_Lambda_NEB_2026_2027_bp:      "Reads classified as Lambda phage 2027bp spike-in",
            fraction_Lambda_NEB_2026_2027_bp: "Fraction of reads classified as Lambda phage 2027bp spike-in",
            num_Lambda_NEB_2026_2322_bp:      "Reads classified as Lambda phage 2322bp spike-in",
            fraction_Lambda_NEB_2026_2322_bp: "Fraction of reads classified as Lambda phage 2322bp spike-in",
            num_Lambda_NEB_2026_4361_bp:      "Reads classified as Lambda phage 4361bp spike-in",
            fraction_Lambda_NEB_2026_4361_bp: "Fraction of reads classified as Lambda phage 4361bp spike-in",
            num_Lambda_NEB_2026_6557_bp:      "Reads classified as Lambda phage 6557bp spike-in",
            fraction_Lambda_NEB_2026_6557_bp: "Fraction of reads classified as Lambda phage 6557bp spike-in",
            num_Lambda_NEB_2026_9416_bp:      "Reads classified as Lambda phage 9416bp spike-in",
            fraction_Lambda_NEB_2026_9416_bp: "Fraction of reads classified as Lambda phage 9416bp spike-in",
            num_Lambda_NEB_2026_23130_bp:     "Reads classified as Lambda phage 23130bp spike-in",
            fraction_Lambda_NEB_2026_23130_bp:"Fraction of reads classified as Lambda phage 23130bp spike-in",
            num_pBR322_BamHI:                 "Reads classified as pBR322 spike-in",
            fraction_pBR322_BamHI:            "Fraction of reads classified as pBR322 spike-in",
            num_phiX174_NEB_PstI:             "Reads classified as phiX174 spike-in",
            fraction_phiX174_NEB_PstI:        "Fraction of reads classified as phiX174 spike-in",
            num_M13mp18_PstI:                 "Reads classified as M13mp18 spike-in",
            fraction_M13mp18_PstI:            "Fraction of reads classified as M13mp18 spike-in",
            kraken_report:                "Kraken2 report (per-taxon read counts and percentages)",
            kraken_output:                "Per-read kraken2 classification",
            kraken2_stats:                "Parsed kraken2 summary TSV (key taxa counts and percentages)",
            pct_bacteria:                 "Percentage of reads classified as Bacteria",
            pct_virus:                    "Percentage of reads classified as Viruses",
            pct_fungi:                    "Percentage of reads classified as Fungi",
            pct_human:                    "Percentage of reads classified as Homo sapiens",
            pct_unclassified:             "Percentage of reads left unclassified",
            top_genus:                    "Genus with the highest kraken2 classification percentage",
            pct_top_genus:                "Percentage of reads assigned to top_genus",
            top_species:                  "Species with the highest kraken2 classification percentage",
            pct_top_species:              "Percentage of reads assigned to top_species"
        }
    }

    parameter_meta {
        input_bam:          "Raw unmapped PacBio HiFi BAM for the sample"
        sample_name:        "Sample identifier used as output file prefix (e.g. bc2097)"
        spikein_fasta:      "Spike-in reference FASTA (NEB HiFi spike-in panel, 11 sequences)"
        kraken2_db_tgz:     "Kraken2 database as a single compressed archive (.tar.zst); extracted at runtime. Use gs://pathogen-public-dbs/jhu/k2_pluspf_20250714.tar.zst for the PlusPF DB."
        kraken2_confidence: "Kraken2 confidence threshold (default 0.001)"
    }

    input {
        File   input_bam
        String sample_name
        File   spikein_fasta
        File   kraken2_db_tgz
        Float  kraken2_confidence = 0.001
    }

    call Prep.BamToFastqAndStats as t_01_BamToFastqAndStats {
        input:
            input_bam = input_bam
    }

    call QC.HifiSeqkitStats as t_02_HifiSeqkitStats {
        input:
            input_fastq = t_01_BamToFastqAndStats.fastq_gz
    }

    call Clean.SpikeInRemoval as t_03_SpikeInRemoval {
        input:
            input_fastq  = t_01_BamToFastqAndStats.fastq_gz,
            spikein_fasta = spikein_fasta
    }

    call QC.HifiKraken2 as t_04_HifiKraken2 {
        input:
            input_fastq    = t_03_SpikeInRemoval.cleaned_fastq,
            kraken2_db_tgz = kraken2_db_tgz,
            confidence     = kraken2_confidence
    }

    call Clean.FilterCleanBam as t_05_FilterCleanBam {
        input:
            input_bam          = input_bam,
            spikein_read_names = t_03_SpikeInRemoval.spikein_read_names,
            kraken_output      = t_04_HifiKraken2.kraken_output,
            sample_name        = sample_name
    }

    output {
        File  cleaned_bam           = t_05_FilterCleanBam.cleaned_bam
        File  cleaned_bam_bai       = t_05_FilterCleanBam.cleaned_bam_bai
        Int   num_spikein_excluded  = t_05_FilterCleanBam.num_spikein_excluded
        Int   num_human_excluded    = t_05_FilterCleanBam.num_human_excluded
        Int   num_reads_cleaned_bam = t_05_FilterCleanBam.num_reads_cleaned_bam

        Float mean_read_accuracy    = t_01_BamToFastqAndStats.mean_read_accuracy
        Float mean_qual_score       = t_01_BamToFastqAndStats.mean_qual_score
        Int   mean_passes           = t_01_BamToFastqAndStats.mean_passes

        File  seqkit_stats_tsv            = t_02_HifiSeqkitStats.seqkit_stats_tsv
        File  seqkit_fx2tab_tsv           = t_02_HifiSeqkitStats.seqkit_fx2tab_tsv
        Int   num_reads                   = t_02_HifiSeqkitStats.num_reads
        Float bases_in_reads              = t_02_HifiSeqkitStats.bases_in_reads
        Int   max_read_length             = t_02_HifiSeqkitStats.max_read_length
        Int   q1_read_length              = t_02_HifiSeqkitStats.q1_read_length
        Int   median_read_length          = t_02_HifiSeqkitStats.median_read_length
        Int   q3_read_length              = t_02_HifiSeqkitStats.q3_read_length
        Int   n50_read_length             = t_02_HifiSeqkitStats.n50_read_length
        Float pct_q20_bases               = t_02_HifiSeqkitStats.pct_q20_bases
        Float pct_q30_bases               = t_02_HifiSeqkitStats.pct_q30_bases
        Float mean_read_gc                = t_02_HifiSeqkitStats.mean_read_gc
        Float bases_in_reads_over_10kb    = t_02_HifiSeqkitStats.bases_in_reads_over_10kb
        Float bases_in_reads_over_20kb    = t_02_HifiSeqkitStats.bases_in_reads_over_20kb

        File  spikein_report        = t_03_SpikeInRemoval.spikein_report
        File  spikein_stats_tsv     = t_03_SpikeInRemoval.spikein_stats_tsv
        Int   total_spikein_reads   = t_03_SpikeInRemoval.total_spikein_reads
        Float fraction_spikein      = t_03_SpikeInRemoval.fraction_spikein

        Int   num_Lambda_NEB_2026_125_bp       = t_03_SpikeInRemoval.num_Lambda_NEB_2026_125_bp
        Float fraction_Lambda_NEB_2026_125_bp  = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_125_bp
        Int   num_Lambda_NEB_2026_564_bp       = t_03_SpikeInRemoval.num_Lambda_NEB_2026_564_bp
        Float fraction_Lambda_NEB_2026_564_bp  = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_564_bp
        Int   num_Lambda_NEB_2026_2027_bp      = t_03_SpikeInRemoval.num_Lambda_NEB_2026_2027_bp
        Float fraction_Lambda_NEB_2026_2027_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_2027_bp
        Int   num_Lambda_NEB_2026_2322_bp      = t_03_SpikeInRemoval.num_Lambda_NEB_2026_2322_bp
        Float fraction_Lambda_NEB_2026_2322_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_2322_bp
        Int   num_Lambda_NEB_2026_4361_bp      = t_03_SpikeInRemoval.num_Lambda_NEB_2026_4361_bp
        Float fraction_Lambda_NEB_2026_4361_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_4361_bp
        Int   num_Lambda_NEB_2026_6557_bp      = t_03_SpikeInRemoval.num_Lambda_NEB_2026_6557_bp
        Float fraction_Lambda_NEB_2026_6557_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_6557_bp
        Int   num_Lambda_NEB_2026_9416_bp      = t_03_SpikeInRemoval.num_Lambda_NEB_2026_9416_bp
        Float fraction_Lambda_NEB_2026_9416_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_9416_bp
        Int   num_Lambda_NEB_2026_23130_bp     = t_03_SpikeInRemoval.num_Lambda_NEB_2026_23130_bp
        Float fraction_Lambda_NEB_2026_23130_bp = t_03_SpikeInRemoval.fraction_Lambda_NEB_2026_23130_bp
        Int   num_pBR322_BamHI                 = t_03_SpikeInRemoval.num_pBR322_BamHI
        Float fraction_pBR322_BamHI            = t_03_SpikeInRemoval.fraction_pBR322_BamHI
        Int   num_phiX174_NEB_PstI             = t_03_SpikeInRemoval.num_phiX174_NEB_PstI
        Float fraction_phiX174_NEB_PstI        = t_03_SpikeInRemoval.fraction_phiX174_NEB_PstI
        Int   num_M13mp18_PstI                 = t_03_SpikeInRemoval.num_M13mp18_PstI
        Float fraction_M13mp18_PstI            = t_03_SpikeInRemoval.fraction_M13mp18_PstI

        File   kraken_report    = t_04_HifiKraken2.kraken_report
        File   kraken_output    = t_04_HifiKraken2.kraken_output
        File   kraken2_stats    = t_04_HifiKraken2.kraken2_stats
        Float  pct_bacteria     = t_04_HifiKraken2.pct_bacteria
        Float  pct_virus        = t_04_HifiKraken2.pct_virus
        Float  pct_fungi        = t_04_HifiKraken2.pct_fungi
        Float  pct_human        = t_04_HifiKraken2.pct_human
        Float  pct_unclassified = t_04_HifiKraken2.pct_unclassified
        String top_genus        = t_04_HifiKraken2.top_genus
        Float  pct_top_genus    = t_04_HifiKraken2.pct_top_genus
        String top_species      = t_04_HifiKraken2.top_species
        Float  pct_top_species  = t_04_HifiKraken2.pct_top_species
    }
}
