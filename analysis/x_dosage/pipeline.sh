#!/usr/bin/env bash
# Host-read sex/karyotype pipeline glue for one sample.
# fastq -> aligned BAM -> mosdepth windows -> rx_sex.py
# Run short-read and long-read paths SEPARATELY (also serves as a cross-platform swap check).
set -euo pipefail

REF=chm13v2.0_maskedY.rCRS.fa      # UCSC-named T2T-CHM13v2.0, Y-PAR masked. NOT the RefSeq GCF_* fasta.
T=${THREADS:-8}
SAMP=${1:?usage: pipeline.sh SAMPLE_ID}

# --- index once ---
# bwa-mem2 index "$REF"
# samtools faidx "$REF"

# --- short reads (Illumina) ---
# Viral reads that don't match human simply stay unmapped; no host/virus pre-split needed.
bwa-mem2 mem -t "$T" "$REF" "${SAMP}_R1.fq.gz" "${SAMP}_R2.fq.gz" \
  | samtools sort -@4 -o "${SAMP}.sr.bam" -
samtools index "${SAMP}.sr.bam"

# --- long reads (ONT: map-ont | PacBio HiFi: map-hifi | CLR: map-pb) ---
minimap2 -t "$T" -ax map-ont "$REF" "${SAMP}.lr.fq.gz" \
  | samtools sort -@4 -o "${SAMP}.lr.bam" -
samtools index "${SAMP}.lr.bam"

# --- windowed depth (MAPQ>=20, 1Mb) -> *.regions.bed.gz ---
mosdepth -t4 -n --by 1000000 -Q 20 "${SAMP}.sr" "${SAMP}.sr.bam"
mosdepth -t4 -n --by 1000000 -Q 20 "${SAMP}.lr" "${SAMP}.lr.bam"

# --- classify (use the platform-matched calibrated constants once TODO #2 lands) ---
python3 rx_sex.py "${SAMP}.sr.regions.bed.gz"
python3 rx_sex.py "${SAMP}.lr.regions.bed.gz"
# Cross-platform swap check: the SR and LR karyotype calls for a sample must agree.

# ---------------------------------------------------------------------------
# Optional QC (not required by the classifier)
#
# 10s smoke test (count-based Rx; self-adjusts to reference lengths, no PAR exclusion):
#   samtools idxstats "${SAMP}.sr.bam" | awk '$1~/^chr([0-9]+|X)$/{
#     n=$3; if($1=="chrX")x=n/$2; else{a+=n; L+=$2}} END{printf "Rx~%.3f\n", x/(a/L)}'
#
# X-heterozygosity corroboration (depth-independent; T2T non-PAR X coords):
#   bcftools mpileup -r chrX:2394411-153925833 -f "$REF" "${SAMP}.sr.bam" \
#     | bcftools call -mv | bcftools view -e 'INFO/DP<8' \
#     | bcftools query -f '[%GT\n]' | sort | uniq -c   # females: clear 0/1 fraction; males ~0
# ---------------------------------------------------------------------------
