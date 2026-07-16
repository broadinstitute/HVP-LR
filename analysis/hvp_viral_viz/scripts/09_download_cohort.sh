#!/usr/bin/env bash
# Run on host with gcloud auth (your identity, not compute SA).
# Pulls 48 foldseek m8 outputs from Terra submission 96f2f736 to local disk.
# (Initially 45 Succeeded on 2026-06-22; final 3 — _4P, _17P, _22P — added
# 2026-06-25 after the submission completed.)
#
# Prereqs:
#   gcloud auth application-default login   (one-time)
#
# Usage:
#   bash scripts/download_cohort.sh
#
# Destination paths below are relative to the analysis/hvp_viral_viz/ project
# root (this script's parent's parent directory).

set -euo pipefail

cd "$(dirname "$0")/.."

DEST=data/cohort_2026-06-22/raw
mkdir -p "$DEST"

gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/9ec0932c-62ea-41c3-ac2f-292b9928e0e4/call-t_08_Format/HVP-0006.1_1P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_1P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/e4bd2497-f5b3-40a4-aacb-422489a8bb59/call-t_08_Format/HVP-0006.1_2P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_2P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/f50f8614-02fb-409e-b5a2-e825d2a43997/call-t_08_Format/HVP-0006.1_3P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_3P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/9c1174fb-78af-4c31-8505-e6034bbb4e05/call-t_08_Format/HVP-0006.1_4P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_4P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/0fb94758-632f-4e03-8ed7-f90d591ed396/call-t_08_Format/HVP-0006.1_5P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_5P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/ef5c1776-17d6-4b31-9722-cdc127a08a75/call-t_08_Format/HVP-0006.1_6P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_6P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/b332fa15-fbf7-4e84-a92a-9cb86839f191/call-t_08_Format/HVP-0006.1_7P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_7P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/0aebd8fc-e0ef-446a-8484-fc74b9a8bed5/call-t_08_Format/HVP-0006.1_8P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_8P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/189c00f5-7744-4e16-a43b-cda539a51fb4/call-t_08_Format/HVP-0006.1_9P.vs_bfvd.m8  "${DEST}"/HVP-0006.1_9P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/78ac5857-5ef6-4fde-99e7-3ad1b9e60dd4/call-t_08_Format/HVP-0006.1_10P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_10P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/75ddd555-fda7-4281-a35d-3aedc6fed140/call-t_08_Format/HVP-0006.1_11P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_11P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/99043861-4208-4289-9b8d-9d4bad621ce1/call-t_08_Format/HVP-0006.1_12P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_12P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/6c21ac3a-73ba-43b3-b6a1-ac98a63f28bb/call-t_08_Format/HVP-0006.1_13P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_13P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/9876040b-db80-4c7b-b248-4da07e61c555/call-t_08_Format/HVP-0006.1_14P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_14P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/736ce7a5-d526-4ada-a57d-ac67a5a9a3a7/call-t_08_Format/HVP-0006.1_15P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_15P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/bec86907-20ef-4871-b157-b24bdd47f65d/call-t_08_Format/HVP-0006.1_16P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_16P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/e1332a54-8e9e-4eb4-892a-548d6f61abc9/call-t_08_Format/HVP-0006.1_17P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_17P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/20373acb-c607-4356-a5bc-3d0eed6ea251/call-t_08_Format/HVP-0006.1_18P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_18P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/2694c304-9729-4313-b67b-ab9b64423457/call-t_08_Format/HVP-0006.1_19P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_19P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/c63e177d-2776-40a1-878c-e055d29214b7/call-t_08_Format/HVP-0006.1_20P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_20P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/24879ee9-1a7f-4d9b-9f48-f8383bc0b113/call-t_08_Format/HVP-0006.1_21P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_21P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/cfefddad-8fb3-4f7b-be00-c6548b15dfc1/call-t_08_Format/HVP-0006.1_22P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_22P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/8233cc1e-4b9a-4722-85e2-df42074b728c/call-t_08_Format/HVP-0006.1_23P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_23P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/ecb3e932-0b9b-44d8-a1bc-d5ed71ecc86c/call-t_08_Format/HVP-0006.1_24P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_24P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/56736622-9ae0-4ece-90f6-229b5d9b55f4/call-t_08_Format/HVP-0006.1_25P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_25P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/d45fced1-514f-4c29-a0bb-3c01880897c6/call-t_08_Format/HVP-0006.1_26P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_26P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/fe586d4d-dfec-4f9d-90c0-dbd72680d4ac/call-t_08_Format/HVP-0006.1_27P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_27P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/fe9d6dd8-3128-40c1-8dfa-c4a65c46abed/call-t_08_Format/HVP-0006.1_28P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_28P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/9e69299e-3643-4ba5-a1f6-e4dbca8c033c/call-t_08_Format/HVP-0006.1_29P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_29P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/c40014be-1ce3-4fa0-a523-c53a052b3499/call-t_08_Format/HVP-0006.1_30P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_30P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/8b581b65-ec49-4a86-85c3-9f325af3d83d/call-t_08_Format/HVP-0006.1_31P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_31P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/4cfb664d-da31-4f53-86d1-641fc222e660/call-t_08_Format/HVP-0006.1_32P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_32P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/100b5b1d-5a37-49ac-838f-a9fe672f9295/call-t_08_Format/HVP-0006.1_33P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_33P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/629d635a-6bc6-489d-849e-d49485c59cfb/call-t_08_Format/HVP-0006.1_34P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_34P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/e0ec5d0d-142a-4831-bae1-4c041bdaa7f2/call-t_08_Format/HVP-0006.1_35P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_35P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/6e13cfb7-b0d0-4d41-b713-87014821c0ea/call-t_08_Format/attempt-2/HVP-0006.1_36P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_36P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/2781d341-30ea-4ddd-bc0e-e6e43d77b911/call-t_08_Format/HVP-0006.1_37P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_37P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/80b4ca64-ae3d-4aa3-b46f-b124c5f2c1a1/call-t_08_Format/HVP-0006.1_38P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_38P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/b40e8e08-45fe-4eb4-8da9-69786147184f/call-t_08_Format/HVP-0006.1_39P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_39P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/2c6e2b63-ffe6-4a51-a70f-eeabe23eaf10/call-t_08_Format/HVP-0006.1_40P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_40P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/c9be4688-baaa-447a-897c-33fd3eae47a2/call-t_08_Format/HVP-0006.1_41P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_41P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/ec042ffe-cc34-4a0d-96cd-3668167104f8/call-t_08_Format/HVP-0006.1_42P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_42P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/8ca76695-fa78-491e-a591-fdc315d3e104/call-t_08_Format/HVP-0006.1_43P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_43P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/88ecd9c1-1460-4ceb-a177-7524c81c7a1b/call-t_08_Format/HVP-0006.1_44P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_44P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/6e3538f6-91ef-46f9-9fd9-5bb86520574c/call-t_08_Format/HVP-0006.1_45P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_45P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/c05801c3-6cd8-4164-93b1-79e310305e2f/call-t_08_Format/HVP-0006.1_46P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_46P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/d004112f-4826-4034-9d64-676baa711be9/call-t_08_Format/HVP-0006.1_47P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_47P.vs_bfvd.m8
gcloud storage cp --no-clobber gs://fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b/submissions/96f2f736-f477-4a62-baab-d43158b9278b/HvpViralProteinAnnotation/ba7f5593-d305-441c-b619-a0cb30acc388/call-t_08_Format/HVP-0006.1_48P.vs_bfvd.m8 "${DEST}"/HVP-0006.1_48P.vs_bfvd.m8
