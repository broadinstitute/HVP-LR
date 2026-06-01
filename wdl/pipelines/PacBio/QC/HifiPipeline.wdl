version 1.0

import "../../../tasks/Utility/Taxonomy.wdl" as Tax
import "../../../tasks/Preprocessing/BamConversion.wdl" as Prep
import "../../../tasks/QC/HifiQC.wdl" as QC

workflow HifiPipeline {

    meta {
        description: "Single-sample PacBio HiFi read QC pipeline. Takes one BAM plus genus/species/sample_name, performs taxonomy lookup, converts BAM to FASTQ, computes seqkit and kraken2 stats, and aggregates per-sample metrics into TSV/TXT/report files."

        allowNestedInputs: true

        outputs: {
            sample_stats_tsv:     "Per-sample summary stats TSV",
            sample_stats_txt:     "Per-sample summary stats in plain-text key:value form",
            sample_stats_report:  "Per-sample human-readable formatted report",
            sample_fastq:         "Gzipped FASTQ produced from the input BAM",
            sample_kraken_report: "Raw kraken2 report",
            sample_taxid_tsv:     "One-row TSV with tax_id, expected_genome_size, genus, species",
            sample_bam_stats:     "Mean read accuracy, Phred quality score, and CCS pass count parsed from BAM tags"
        }
    }

    parameter_meta {
        input_bam:          "PacBio HiFi reads BAM file for the sample"
        sample_name:        "Sample identifier used as the output prefix for per-sample artifacts (e.g. bc2019)"
        genus:              "Genus name. Empty string allowed for metagenomic samples (but disables taxonomy/coverage estimation)."
        species:            "Species name. Empty string allowed for metagenomic samples."
        kraken2_db_hash:    "Kraken2 database hash.k2d file (pre-extracted)"
        kraken2_db_opts:    "Kraken2 database opts.k2d file (pre-extracted)"
        kraken2_db_taxo:    "Kraken2 database taxo.k2d file (pre-extracted)"
        kraken2_confidence: "Kraken2 confidence threshold"
    }

    input {
        File   input_bam
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

    call QC.HifiReadStats as t_05_HifiReadStats {
        input:
            sample_name               = sample_name,
            bam_stats_tsv             = t_02_BamToFastqAndStats.bam_stats_tsv,
            seqkit_stats_tsv          = t_03_HifiSeqkitStats.seqkit_stats_tsv,
            seqkit_fx2tab_tsv         = t_03_HifiSeqkitStats.seqkit_fx2tab_tsv,
            kraken2_stats_tsv         = t_04_HifiKraken2.kraken2_stats,
            taxid_and_genome_size_tsv = t_01_GetTaxIdAndGenomeSize.taxid_and_genome_size_tsv
    }

    output {
        File sample_stats_tsv     = t_05_HifiReadStats.hifi_read_stats_tsv
        File sample_stats_txt     = t_05_HifiReadStats.hifi_read_stats_txt
        File sample_stats_report  = t_05_HifiReadStats.hifi_read_stats_report
        File sample_fastq         = t_02_BamToFastqAndStats.fastq_gz
        File sample_kraken_report = t_04_HifiKraken2.kraken_report
        File sample_taxid_tsv     = t_01_GetTaxIdAndGenomeSize.taxid_and_genome_size_tsv
        File sample_bam_stats     = t_02_BamToFastqAndStats.bam_stats_tsv
    }
}
