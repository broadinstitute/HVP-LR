version 1.0

import "../../../tasks/Preprocessing/BamConversion.wdl" as Prep
import "../../../tasks/Assembly/Myloasm.wdl"             as Asm
import "../../../tasks/Alignment/Minimap2AlignReads.wdl" as Aln

workflow HvpAssembly {

    meta {
        description: "Single-sample HVP metagenomic assembly workflow. Accepts a cleaned (spike-in-removed, human-depleted) BAM, converts it to FASTQ internally, assembles with myloasm, aligns reads back to contigs with minimap2, and computes per-contig assembly statistics with seqkit. BAM format is preserved as the upstream currency; FASTQ conversion happens inside this workflow and is not exposed as an output. Returns the primary assembly FASTA, read-to-contig BAM, and scalar QC metrics. The read-to-contig BAM produced here is consumed by the downstream HvpReadRescue workflow."

        allowNestedInputs: true

        outputs: {
            assembly_primary_fa:  "Primary assembly FASTA from myloasm",
            sorted_bam:           "Read-to-contig alignment BAM (minimap2 map-hifi, coordinate-sorted)",
            sorted_bam_bai:       "BAI index for sorted_bam",
            flagstat_txt:         "samtools flagstat output for the read-to-contig alignment",
            align_stats_tsv:      "Alignment summary TSV (total_reads, mapped_reads, fraction_mapped)",
            mapped_reads:         "Total primary reads mapped to assembly contigs",
            fraction_mapped:      "Fraction of primary reads mapped to assembly contigs",
            mean_read_accuracy:   "Mean read accuracy from rq:f BAM tags (percentage)",
            mean_qual_score:      "Mean Phred quality score derived from mean_read_accuracy",
            mean_passes:          "Mean number of CCS subreads per read from np:i BAM tags",
            asm_stats_tsv:        "Assembly contig statistics TSV (12 metrics, header row + values row)",
            num_contigs:          "Number of contigs in the assembly",
            bases_in_contigs:     "Total bases across all contigs",
            mean_contig_length:   "Mean contig length",
            q1_contig_length:     "Q1 (25th percentile) contig length",
            median_contig_length: "Median contig length",
            q3_contig_length:     "Q3 (75th percentile) contig length",
            n50_contig_length:    "N50 contig length",
            max_contig_length:    "Maximum contig length",
            mean_contig_gc:       "Mean GC content (percentage) across contigs",
            num_circ_contigs:     "Number of circular contigs",
            num_1Mb_contigs:      "Number of contigs longer than 1 Mb",
            num_circ_1Mb_contigs: "Number of contigs that are both circular and longer than 1 Mb"
        }
    }

    parameter_meta {
        input_bam:      "Cleaned PacBio HiFi BAM (spike-in-removed, human-depleted) produced by HvpReadProcessing"
        sample_name:    "Sample identifier used as output file prefix (e.g. bc2097)"
        min_contig_len: "Minimum contig length for seqkit assembly stats filtering; 0 = no filtering"
    }

    input {
        File   input_bam
        String sample_name
        Int    min_contig_len = 0
    }

    call Prep.BamToFastqAndStats as t_01_BamToFastqAndStats {
        input:
            input_bam = input_bam
    }

    call Asm.Myloasm as t_02_Myloasm {
        input:
            input_fastq = t_01_BamToFastqAndStats.fastq_gz
    }

    call Aln.Minimap2AlignReads as t_03_Minimap2AlignReads {
        input:
            input_fastq    = t_01_BamToFastqAndStats.fastq_gz,
            assembly_fasta = t_02_Myloasm.assembly_primary_fa,
            sample_name    = sample_name
    }

    call Asm.SeqkitAssemblyStats as t_04_SeqkitAssemblyStats {
        input:
            assembly_fasta = t_02_Myloasm.assembly_primary_fa,
            min_contig_len = min_contig_len
    }

    output {
        File  assembly_primary_fa  = t_02_Myloasm.assembly_primary_fa

        File  sorted_bam           = t_03_Minimap2AlignReads.sorted_bam
        File  sorted_bam_bai       = t_03_Minimap2AlignReads.sorted_bam_bai
        File  flagstat_txt         = t_03_Minimap2AlignReads.flagstat_txt
        File  align_stats_tsv      = t_03_Minimap2AlignReads.align_stats_tsv
        Int   mapped_reads         = t_03_Minimap2AlignReads.mapped_reads
        Float fraction_mapped      = t_03_Minimap2AlignReads.fraction_mapped

        Float mean_read_accuracy   = t_01_BamToFastqAndStats.mean_read_accuracy
        Float mean_qual_score      = t_01_BamToFastqAndStats.mean_qual_score
        Int   mean_passes          = t_01_BamToFastqAndStats.mean_passes

        File  asm_stats_tsv        = t_04_SeqkitAssemblyStats.asm_stats_tsv
        Int   num_contigs          = t_04_SeqkitAssemblyStats.num_contigs
        Int   bases_in_contigs     = t_04_SeqkitAssemblyStats.bases_in_contigs
        Float mean_contig_length   = t_04_SeqkitAssemblyStats.mean_contig_length
        Int   q1_contig_length     = t_04_SeqkitAssemblyStats.q1_contig_length
        Int   median_contig_length = t_04_SeqkitAssemblyStats.median_contig_length
        Int   q3_contig_length     = t_04_SeqkitAssemblyStats.q3_contig_length
        Int   n50_contig_length    = t_04_SeqkitAssemblyStats.n50_contig_length
        Int   max_contig_length    = t_04_SeqkitAssemblyStats.max_contig_length
        Float mean_contig_gc       = t_04_SeqkitAssemblyStats.mean_contig_gc
        Int   num_circ_contigs     = t_04_SeqkitAssemblyStats.num_circ_contigs
        Int   num_1Mb_contigs      = t_04_SeqkitAssemblyStats.num_1Mb_contigs
        Int   num_circ_1Mb_contigs = t_04_SeqkitAssemblyStats.num_circ_1Mb_contigs
    }
}
