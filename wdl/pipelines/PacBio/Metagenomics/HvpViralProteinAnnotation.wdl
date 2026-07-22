version 1.0

import "../../../tasks/ProteinAnnotation/PyrodigalGv.wdl"                  as PG
import "../../../tasks/ProteinAnnotation/Mmseqs2.wdl"                      as MM
import "../../../tasks/ProteinAnnotation/Foldseek.wdl"                     as FS
import "../../../tasks/ProteinAnnotation/ProteinAnnotationHelpers.wdl"     as PAH

workflow HvpViralProteinAnnotation {

    meta {
        description: "Per-sample structural protein annotation for HVP viral contigs. Calls ORFs on every nucleotide source (VS2 viral combined FASTA, optional assembly contigs, optional rescued reads) with Pyrodigal-gv, concatenates the predicted proteins with geNomad's pre-computed virus_proteins.faa under per-source header prefixes, collapses the union to a non-redundant representative set with MMseqs2 easy-linclust (default 90% AA identity, 80% coverage), folds the NR proteins into a foldseek structure DB via ProstT5, structurally searches that DB against BFVD with foldseek search, formats the alignment as a TSV with foldseek convertalis, and (optionally) joins each query's best hit against a BFVD reference-metadata TSV to attach descriptive annotations. Inputs that are gzipped are auto-decompressed by each tool. All major intermediates are surfaced as separate File outputs so downstream tooling can pick its level."

        allowNestedInputs: true

        outputs: {
            vs2_proteins_faa:         "Pyrodigal-gv predicted proteins (amino-acid FASTA) called on the VS2 viral_combined_fa nucleotide input.",
            assembly_proteins_faa:    "Pyrodigal-gv predicted proteins called on assembly_contigs_fa, when supplied. Absent if no assembly_contigs_fa input.",
            rescued_proteins_faa:     "Pyrodigal-gv predicted proteins called on rescued_reads_fa_gz, when supplied. Absent if no rescued_reads_fa_gz input.",
            combined_proteins_faa:    "All input AA FASTAs concatenated with per-source header prefixes (vs2|, genomad|, assembly|, rescued|).",
            nr_proteins_faa:          "Non-redundant representative protein set after mmseqs easy-linclust (one sequence per cluster).",
            nr_cluster_tsv:           "Two-column cluster-membership TSV from mmseqs (representative <TAB> member) for the NR collapse.",
            foldseek_db_archive:      "tar.gz of the foldseek structure DB built from nr_proteins_faa via ProstT5 (3Di tokens predicted from AA).",
            foldseek_aln_db_archive:  "tar.gz of the foldseek search alignment DB against BFVD.",
            foldseek_hits_tsv:        "BLAST-tab-like TSV of all foldseek hits (columns: query,target,evalue,bits,fident,alnlen,prob).",
            best_hits_tsv:            "One row per query — lowest-e-value foldseek hit, ties broken by higher bit score.",
            annotated_hits_tsv:       "best_hits_tsv left-joined on target against bfvd_metadata_tsv (when supplied). If no metadata supplied, identical schema to best_hits_tsv.",
            num_vs2_orfs:             "ORFs called on VS2 viral combined FASTA.",
            num_assembly_orfs:        "ORFs called on assembly contigs (0 if no assembly input supplied).",
            num_rescued_orfs:         "ORFs called on rescued reads (0 if no rescued-reads input supplied).",
            num_combined_proteins:    "Total sequences in combined_proteins_faa before NR collapse.",
            num_nr_proteins:          "Representative protein count after mmseqs easy-linclust.",
            num_foldseek_hits:        "Total foldseek hit rows in foldseek_hits_tsv (query-target pairs).",
            num_queries_with_hits:    "Distinct query proteins with at least one foldseek hit.",
            num_queries_annotated:    "Best-hit queries whose target matched a row in bfvd_metadata_tsv. Zero if no metadata supplied."
        }
    }

    parameter_meta {
        sample_name:               "Sample identifier used as output file prefix (e.g. bc2097)."
        vs2_viral_combined_fa:     "VirSorter2 viral_combined_fa nucleotide FASTA (output of the VS2 task in HvpViralPipeline)."
        genomad_virus_proteins_faa: "geNomad virus_proteins.faa amino-acid FASTA. Concatenated directly — no ORF call needed."
        assembly_contigs_fa:       "Optional assembly contigs FASTA (e.g. metaFlye contigs) to ORF-call separately for novel-virus discovery in non-VS2/geNomad-called material."
        rescued_reads_fa_gz:       "Optional rescued-reads FASTA (gzipped, .fa.gz), e.g. HvpReadRescue rescue_viral_fa_gz, to ORF-call separately."
        bfvd_db_tgz:               "Pre-built BFVD foldseek structure DB packaged as tar.gz. Source: https://bfvd.steineggerlab.workers.dev/latest/bfvd_foldseekdb.tar.gz (~534 MB; stage once in Terra workspace data; reuse across runs). Version manifest: https://bfvd.steineggerlab.workers.dev/latest/bfvd.version"
        bfvd_db_tgz_padded:        "Optional padded variant of bfvd_db_tgz, required for use_gpu=true. Produced once on a GPU workstation via foldseek 10 makepaddedseqdb (see FoldseekSearch comment block for the recipe; ~14 GB output, MUST be roughly the same size as the unpadded source). When use_gpu=true and this is supplied, FoldseekSearch consumes the padded tarball; when use_gpu=false, this input is ignored and the unpadded bfvd_db_tgz is used as-is."
        prostt5_weights_tgz:       "ProstT5 model weights packaged as tar.gz. Source: https://foldseek.steineggerlab.workers.dev/prostt5-f16-gguf.tar.gz (~2.1 GB)."
        bfvd_metadata_tsv:         "Optional BFVD reference metadata TSV (tab-delimited, headerless, 6 columns: UniRef100/UniProt accession, model filename, avg_pLDDT, pTM, splitted flag, version). First column joins against foldseek target IDs. Source: https://bfvd.steineggerlab.workers.dev/latest/bfvd_metadata.tsv (~15 MB). Leave unset to skip annotation join. Multiple rows per UniProt key (one per structure-model split) collapse to last-row-wins in the join."
        mmseqs_min_seq_id:         "Minimum amino-acid identity for the mmseqs easy-linclust NR collapse. Default 0.9."
        mmseqs_coverage:           "Minimum alignment coverage for the mmseqs easy-linclust NR collapse. Default 0.8."
        foldseek_evalue:           "Maximum e-value for foldseek search hits. Default 0.001."
        pyrodigal_extra_args:      "Additional command-line args appended verbatim to every Pyrodigal-gv invocation."
        foldseek_extra_args:       "Additional command-line args appended verbatim to the foldseek search invocation."
        use_gpu:                   "When true (default), run FoldseekCreateDbFromFasta + FoldseekSearch on GPU (foldseek-gpu image, nvidia-tesla-t4). Requires bfvd_db_tgz_padded to be supplied (FoldseekSearch will fail at TPRE detection otherwise). When false, both tasks run CPU-only on the foldseek image and the unpadded bfvd_db_tgz is used directly."
    }

    input {
        String sample_name

        File   vs2_viral_combined_fa
        File   genomad_virus_proteins_faa

        File?  assembly_contigs_fa
        File?  rescued_reads_fa_gz

        File   bfvd_db_tgz
        File?  bfvd_db_tgz_padded
        File   prostt5_weights_tgz
        File?  bfvd_metadata_tsv

        Float   mmseqs_min_seq_id    = 0.9
        Float   mmseqs_coverage      = 0.8
        Float   foldseek_evalue      = 0.001
        String  pyrodigal_extra_args = ""
        String  foldseek_extra_args  = ""
        Boolean use_gpu              = true
    }

    # 1. ORF-call every nucleotide source. VS2 always; assembly + rescued reads only if supplied.
    call PG.PyrodigalGvCallOrfs as t_01_OrfsVS2 {
        input:
            nucleotide_fasta = vs2_viral_combined_fa,
            sample_name      = sample_name,
            source_label     = "vs2",
            extra_args       = pyrodigal_extra_args
    }

    if (defined(assembly_contigs_fa)) {
        call PG.PyrodigalGvCallOrfs as t_02_OrfsAssembly {
            input:
                nucleotide_fasta = select_first([assembly_contigs_fa]),
                sample_name      = sample_name,
                source_label     = "assembly",
                extra_args       = pyrodigal_extra_args
        }
    }

    if (defined(rescued_reads_fa_gz)) {
        call PG.PyrodigalGvCallOrfs as t_03_OrfsRescued {
            input:
                nucleotide_fasta = select_first([rescued_reads_fa_gz]),
                sample_name      = sample_name,
                source_label     = "rescued",
                extra_args       = pyrodigal_extra_args
        }
    }

    # 2. Concatenate every AA FASTA (geNomad pre-computed + per-source Pyrodigal-gv).
    #    The helper task tags every header with its source label so post-cluster hits remain attributable.
    Array[File]   all_fastas = select_all([
        genomad_virus_proteins_faa,
        t_01_OrfsVS2.proteins_faa,
        t_02_OrfsAssembly.proteins_faa,
        t_03_OrfsRescued.proteins_faa
    ])
    Array[String] all_labels = if defined(assembly_contigs_fa) then
                                  (if defined(rescued_reads_fa_gz)
                                   then ["genomad", "vs2", "assembly", "rescued"]
                                   else ["genomad", "vs2", "assembly"])
                              else
                                  (if defined(rescued_reads_fa_gz)
                                   then ["genomad", "vs2", "rescued"]
                                   else ["genomad", "vs2"])

    call PAH.ConcatProteinFastas as t_04_Concat {
        input:
            protein_fastas = all_fastas,
            source_labels  = all_labels,
            prefix         = sample_name
    }

    # 3. Collapse to non-redundant representatives via mmseqs easy-linclust.
    call MM.MmseqsEasyLinclust as t_05_NR {
        input:
            input_fasta = t_04_Concat.combined_faa,
            prefix      = sample_name + ".nr",
            min_seq_id  = mmseqs_min_seq_id,
            coverage    = mmseqs_coverage
    }

    # 4. Fold NR proteins into a foldseek structure DB via ProstT5.
    call FS.FoldseekCreateDbFromFasta as t_06_FoldDb {
        input:
            protein_fasta       = t_05_NR.rep_seq_fasta,
            prefix              = sample_name + ".nr",
            prostt5_weights_tgz = prostt5_weights_tgz,
            use_gpu             = use_gpu
    }

    # 5. Structural search against BFVD.
    #    GPU search requires the padded BFVD DB (foldseek 10 makepaddedseqdb).
    #    When use_gpu=true, prefer bfvd_db_tgz_padded if supplied; if it is NOT
    #    supplied, fall through to the unpadded bfvd_db_tgz — FoldseekSearch's
    #    own TPRE-detect guard will then fail with an explicit error. When
    #    use_gpu=false, the unpadded tarball is the correct input.
    File foldseek_target_db = if use_gpu
        then select_first([bfvd_db_tgz_padded, bfvd_db_tgz])
        else bfvd_db_tgz

    call FS.FoldseekSearch as t_07_Search {
        input:
            query_db_archive  = t_06_FoldDb.db_archive,
            target_db_archive = foldseek_target_db,
            prefix            = sample_name + ".vs_bfvd",
            evalue_cutoff     = foldseek_evalue,
            use_gpu           = use_gpu,
            extra_args        = foldseek_extra_args
    }

    # 6. Format the alignment DB as a TSV.
    #    Note: queries here are ProstT5-built (sequence-only) foldseek DBs — they have no CA
    #    (C-alpha) datafile, so any --format-output column that needs structural coordinates
    #    will fail. Use a CA-free column set; substitute `mismatch` for the structural-similarity
    #    `prob` column to keep a useful 7-column output.
    String hits_format_output_str   = "query,target,evalue,bits,fident,alnlen,mismatch"
    Array[String] hits_columns_list = ["query", "target", "evalue", "bits", "fident", "alnlen", "mismatch"]

    call FS.FoldseekConvertAlis as t_08_Format {
        input:
            query_db_archive  = t_06_FoldDb.db_archive,
            target_db_archive = foldseek_target_db,
            aln_db_archive    = t_07_Search.aln_db_archive,
            format_output     = hits_format_output_str,
            prefix            = sample_name + ".vs_bfvd"
    }

    # 7. Best-hit per query and (optional) annotation transfer from BFVD metadata.
    #    BFVD metadata is headerless; supply explicit column names so the join output
    #    has descriptive headers (target_meta_uniprot_id, target_meta_avg_plddt, ...).
    call PAH.AnnotationTransfer as t_09_Annotate {
        input:
            foldseek_hits_tsv             = t_08_Format.results_tsv,
            hits_columns                  = hits_columns_list,
            reference_metadata_tsv        = bfvd_metadata_tsv,
            reference_metadata_has_header = false,
            reference_metadata_columns    = ["uniprot_id", "model_id", "avg_plddt", "ptm", "splitted", "version"],
            prefix                        = sample_name + ".vs_bfvd"
    }

    output {
        File  vs2_proteins_faa        = t_01_OrfsVS2.proteins_faa
        File? assembly_proteins_faa   = t_02_OrfsAssembly.proteins_faa
        File? rescued_proteins_faa    = t_03_OrfsRescued.proteins_faa

        File  combined_proteins_faa   = t_04_Concat.combined_faa
        File  nr_proteins_faa         = t_05_NR.rep_seq_fasta
        File  nr_cluster_tsv          = t_05_NR.cluster_tsv

        File  foldseek_db_archive     = t_06_FoldDb.db_archive
        File  foldseek_aln_db_archive = t_07_Search.aln_db_archive
        File  foldseek_hits_tsv       = t_08_Format.results_tsv

        File  best_hits_tsv           = t_09_Annotate.best_hits_tsv
        File  annotated_hits_tsv      = t_09_Annotate.annotated_hits_tsv

        Int   num_vs2_orfs            = t_01_OrfsVS2.num_genes
        Int   num_assembly_orfs       = select_first([t_02_OrfsAssembly.num_genes, 0])
        Int   num_rescued_orfs        = select_first([t_03_OrfsRescued.num_genes, 0])
        Int   num_combined_proteins   = t_04_Concat.num_seqs
        Int   num_nr_proteins         = t_05_NR.num_clusters
        Int   num_foldseek_hits       = t_08_Format.num_hits
        Int   num_queries_with_hits   = t_09_Annotate.num_queries_hit
        Int   num_queries_annotated   = t_09_Annotate.num_queries_annotated
    }
}
