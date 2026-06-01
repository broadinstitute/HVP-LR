version 1.0

import "../../../tasks/Utility/SampleSheetUtils.wdl" as Utils
import "../../../tasks/Utility/Taxonomy.wdl" as Tax
import "../../../tasks/Preprocessing/BamConversion.wdl" as Prep
import "../../../tasks/QC/HifiQC.wdl" as QC

workflow HifiPipeline {

    meta {
        description: "PacBio HiFi read QC pipeline. Takes parallel per-sample arrays (BAM path, genus, species, strain, barcode), performs batch taxonomy lookup, then per-sample converts BAM to FASTQ, computes seqkit and kraken2 stats, aggregates metrics per sample, and merges into a single multi-sample report."

        allowNestedInputs: true

        outputs: {
            all_samples_stats_tsv:    "Multi-sample concatenated stats TSV with bam and strain columns appended",
            per_sample_stats_tsv:     "Per-sample summary stats TSV (one per BAM)",
            per_sample_stats_txt:     "Per-sample summary stats in plain-text key:value form",
            per_sample_stats_report:  "Per-sample human-readable formatted report",
            per_sample_fastq:         "Per-sample gzipped FASTQ produced from input BAM",
            per_sample_kraken_report: "Per-sample raw kraken2 report"
        }
    }

    parameter_meta {
        bams:               "Per-sample BAM paths (typically gs:// URIs). One entry per sample; positionally aligned with genera/species/strains/barcodes."
        genera:             "Per-sample genus names. Empty string allowed for metagenomic samples."
        species:            "Per-sample species names. Empty string allowed for metagenomic samples."
        strains:            "Per-sample strain identifiers."
        barcodes:           "Per-sample barcode identifiers (e.g. bc2019), used as the output prefix for per-sample artifacts."
        kraken2_db_hash:    "Kraken2 database hash.k2d file (pre-extracted)"
        kraken2_db_opts:    "Kraken2 database opts.k2d file (pre-extracted)"
        kraken2_db_taxo:    "Kraken2 database taxo.k2d file (pre-extracted)"
        kraken2_confidence: "Kraken2 confidence threshold"
    }

    input {
        Array[String] bams
        Array[String] genera
        Array[String] species
        Array[String] strains
        Array[String] barcodes

        File  kraken2_db_hash    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/hash.k2d"  # !FileCoercion
        File  kraken2_db_opts    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/opts.k2d"  # !FileCoercion
        File  kraken2_db_taxo    = "gs://gcid-cil-shed-archive/kraken_db/k2_pluspf_20251015/taxo.k2d"  # !FileCoercion
        Float kraken2_confidence = 0.001
    }

    call Tax.BatchGetTaxIdAndGenomeSize as t_01_BatchGetTaxIdAndGenomeSize {
        input:
            genera       = genera,
            species_list = species,
            barcodes     = barcodes
    }

    scatter (idx in range(length(bams))) {

        String bam_path = bams[idx]
        String barcode  = barcodes[idx]

        call Prep.BamToFastqAndStats as t_02_BamToFastqAndStats {
            input:
                input_bam = bam_path  # !FileCoercion
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
                sample_name               = barcode,
                bam_stats_tsv             = t_02_BamToFastqAndStats.bam_stats_tsv,
                seqkit_stats_tsv          = t_03_HifiSeqkitStats.seqkit_stats_tsv,
                seqkit_fx2tab_tsv         = t_03_HifiSeqkitStats.seqkit_fx2tab_tsv,
                kraken2_stats_tsv         = t_04_HifiKraken2.kraken2_stats,
                taxid_and_genome_size_tsv = t_01_BatchGetTaxIdAndGenomeSize.taxid_and_genome_size_tsvs[idx]
        }
    }

    call Utils.MergeSampleStats as t_06_MergeSampleStats {
        input:
            sample_tsvs = t_05_HifiReadStats.hifi_read_stats_tsv,
            bam_paths   = bams,
            strains     = strains
    }

    output {
        File        all_samples_stats_tsv    = t_06_MergeSampleStats.merged_stats_tsv
        Array[File] per_sample_stats_tsv     = t_05_HifiReadStats.hifi_read_stats_tsv
        Array[File] per_sample_stats_txt     = t_05_HifiReadStats.hifi_read_stats_txt
        Array[File] per_sample_stats_report  = t_05_HifiReadStats.hifi_read_stats_report
        Array[File] per_sample_fastq         = t_02_BamToFastqAndStats.fastq_gz
        Array[File] per_sample_kraken_report = t_04_HifiKraken2.kraken_report
    }
}
