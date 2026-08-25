version 1.0

import "../../../tasks/QC/Mosdepth.wdl"        as MD
import "../../../tasks/QC/RxSexKaryotype.wdl"  as RX

workflow SexChromKaryotype {

    meta {
        description: "Sex-chromosome karyotype QC for a single sample from host coverage, alignment-reference T2T-CHM13v2.0. Takes a coordinate-sorted BAM (Illumina or PacBio HiFi reads aligned to chm13v2.0_maskedY.rCRS), computes 1 Mb windowed depth with mosdepth, then calls the karyotype from X/Y dosage with the x_dosage rx_sex classifier. Returns the predicted karyotype and a confidence score plus the underlying dosage metrics. Runs identically on short- and long-read data; agreement between a sample's SR and LR calls doubles as a cross-platform swap check. Aneuploidy classes are dosage-based; off-multiple / mosaic samples are reported as OTHER by design."

        outputs: {
            karyotype_call:    "Predicted sex-chromosome karyotype (46,XX / 46,XY / 47,XXY / 47,XYY / 47,XXX / 45,X / 48,XXYY / OTHER / INSUFFICIENT)",
            confidence:        "Posterior probability of the predicted karyotype",
            rx:                "Rx dosage (chrX non-PAR / autosome)",
            ry:                "Ry dosage (chrY euchromatin / autosome)",
            auto_depth:        "Median autosomal window depth",
            runner_up:         "Second-ranked karyotype and its posterior (label:prob)",
            karyotype_tsv:     "Full classifier output TSV",
            regions_bed:       "mosdepth 1 Mb windowed depth BED",
            mosdepth_summary:  "mosdepth per-chromosome summary"
        }
    }

    parameter_meta {
        aligned_bam:  "Coordinate-sorted BAM aligned to T2T-CHM13v2.0 (chm13v2.0_maskedY.rCRS)"
        aligned_bai:  "BAI index for aligned_bam"
        sample_name:  "Sample identifier; used as the output prefix"
        platform:     "Calibrated config to use: short_read | long_read (default short_read)"
        config:       "Optional TOML config override; defaults to the baked per-platform config in the x_dosage image"
        window_size:  "mosdepth window size in bp (default 1000000 — the classifier's calibrated granularity)"
        min_mapq:     "Minimum MAPQ for mosdepth (default 20)"
    }

    input {
        File aligned_bam
        File aligned_bai
        String sample_name

        String platform = "short_read"
        File? config

        Int window_size = 1000000
        Int min_mapq = 20
    }

    call MD.Mosdepth as t_01_Mosdepth {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            prefix      = sample_name,
            window_size = window_size,
            min_mapq    = min_mapq
    }

    call RX.RxSexKaryotype as t_02_RxSexKaryotype {
        input:
            regions_bed = t_01_Mosdepth.regions_bed,
            sample_name = sample_name,
            platform    = platform,
            config      = config
    }

    output {
        String karyotype_call = t_02_RxSexKaryotype.karyotype_call
        Float confidence      = t_02_RxSexKaryotype.confidence
        Float rx              = t_02_RxSexKaryotype.rx
        Float ry              = t_02_RxSexKaryotype.ry
        Float auto_depth      = t_02_RxSexKaryotype.auto_depth
        String runner_up      = t_02_RxSexKaryotype.runner_up
        File karyotype_tsv    = t_02_RxSexKaryotype.karyotype_tsv
        File regions_bed      = t_01_Mosdepth.regions_bed
        File mosdepth_summary = t_01_Mosdepth.summary_txt
    }
}
