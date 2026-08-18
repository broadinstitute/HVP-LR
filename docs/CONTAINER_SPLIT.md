# HVP-LR container split — task → tools inventory

Every WDL task, the tools it actually invokes, and the image it runs on today.
`hvp-monolith:0.0.3` is the fat image being retired. "✓ already small" = task
already points at a focused public or custom image and needs no change.

## Tasks already on small images (no change needed)

| Task | File | Tools | Current image |
|------|------|-------|---------------|
| BamToFastqAndStats | Preprocessing/BamConversion | samtools, awk | quay biocontainers/samtools:1.23 ✓ |
| FilterCleanBam | Preprocessing/ReadCleaning | samtools, awk | quay biocontainers/samtools:1.23 ✓ |
| HifiSeqkitStats | QC/HifiQC | seqkit, awk | staphb/seqkit:2.12.0 ✓ |
| SeqkitAssemblyStats | Assembly/Myloasm | seqkit, awk | staphb/seqkit:2.12.0 ✓ |
| CreateHifiQCReport | QC/HifiQC | bash, awk | ubuntu:22.04 ✓ |
| GetTaxIdAndGenomeSize | Utility/Taxonomy | taxonkit | quay biocontainers/taxonkit:0.20 ✓ |
| PyrodigalGvCallOrfs | ProteinAnnotation/PyrodigalGv | pyrodigal-gv | pyrodigal-gv:0.3.2 ✓ |
| MmseqsEasy* / Mmseqs* (13 tasks) | ProteinAnnotation/Mmseqs2 | mmseqs | mmseqs2:18.0.0 ✓ |
| Foldseek* (18 tasks) | ProteinAnnotation/Foldseek | foldseek | foldseek:10.0.1 / foldseek-gpu:10.0.1 ✓ |

## Tasks on the monolith (must be rehomed)

| Task | File | Tools invoked | Monolith env | New image |
|------|------|---------------|--------------|-----------|
| SpikeInRemoval | Preprocessing/ReadCleaning | minimap2, samtools, seqkit, awk | base | **hvp-align** |
| Minimap2AlignReads | Alignment/Minimap2AlignReads | minimap2, samtools, awk | base | **hvp-align** |
| Myloasm | Assembly/Myloasm | myloasm, seqkit, awk | base | **myloasm** |
| HifiKraken2 | QC/HifiQC | kraken2, awk | base | **kraken2** |
| JgiDepth | Binning/Binning | jgi_summarize_bam_contig_depths, awk | base (metabat2) | **hvp-binning** |
| MetaBAT2 | Binning/Binning | metabat2, awk | base | **hvp-binning** |
| MaxBin2 | Binning/Binning | run_MaxBin.pl (maxbin2), awk | base | **hvp-binning** |
| SemiBin2 | Binning/Binning | SemiBin2, awk | base | **hvp-binning** |
| DASTool | Binning/DASTool | DAS_Tool, awk | base | **hvp-binning** |
| CheckM2 | MAG/CheckM2 | checkm2, diamond, awk | checkm2 (py3.10) | **checkm2** |
| Skani | MAG/Skani | skani, awk | base | **skani** |
| Genomad | Viral/Genomad | genomad, awk | base | **genomad** |
| CheckV | Viral/CheckV | checkv, (prodigal), awk | base | **checkv** |
| VirSorter2 | Viral/VirSorter2 | virsorter, python3 | viral (py3.10) | **virsorter2** |
| SkaniAnnotate | MAG/Skani | python3, awk | base (glue) | **hvp-pyutils** |
| BinSummary | MAG/MagSummary | python3, pandas | base (glue) | **hvp-pyutils** |
| MagSummary | MAG/MagSummary | python3, pandas, seqkit | base (glue) | **hvp-pyutils** |
| ReadRescue | Rescue/ReadRescue | python3, pysam, samtools, seqkit, awk | base (glue) | **hvp-pyutils** |
| MergeContigsRescue | Rescue/MergeContigsRescue | seqkit, awk | base (glue) | **hvp-pyutils** |
| ConcatProteinFastas | ProteinAnnotation/ProteinAnnotationHelpers | python3, seqkit, awk | base (glue) | **hvp-pyutils** |
| AnnotationTransfer | ProteinAnnotation/ProteinAnnotationHelpers | python3, pandas, awk | base (glue) | **hvp-pyutils** |
| ViralContigSummary | Viral/ViralSummary | python3, pandas | base (glue) | **hvp-pyutils** |
| ViralOverallSummary | Viral/ViralSummary | python3, pandas | base (glue) | **hvp-pyutils** |

## Dead weight in the monolith (no WDL task calls these) — drop entirely

`deepvirfinder`, `vcontact2`, `vibrant`, `flye`, `hifiasm_meta`, `metamdbg`,
`concoct`, `bracken`, `fastqc`, `multiqc`, standalone `diamond`.
Removing these is most of the size win and eliminates the py3.6 (deepvirfinder)
and py3.8 (vcontact2) sidecar envs the monolith juggles with wrapper shims.

## Proposed new images (10)

| Image | Tools | Base | Why grouped |
|-------|-------|------|-------------|
| hvp-align   | minimap2 + samtools + seqkit | micromamba | shared by 2 alignment/cleaning tasks; all py3.12-compatible |
| myloasm     | myloasm + seqkit | micromamba | single assembler |
| kraken2     | kraken2 | micromamba/public | single classifier |
| hvp-binning | metabat2 + maxbin2 + semibin2 + das_tool | micromamba | 5 binning tasks, one conda solve, all py3.12 |
| checkm2     | checkm2 + diamond | micromamba (py3.10) | own python pin — cannot share base |
| skani       | skani | micromamba | single ANI tool |
| genomad     | genomad | micromamba | single viral caller, large DB |
| checkv      | checkv | micromamba | single QC tool, own DB |
| virsorter2  | virsorter2 | micromamba (py3.10) | own python pin — cannot share base |
| hvp-pyutils | python3.12 + pandas + numpy + pysam + biopython + seqkit | micromamba | all pure-python glue/summary tasks (9) |

hvp-pyutils absorbs every task whose "tool" is just pandas/pysam scripting, so
those never drag a heavy tool image.
