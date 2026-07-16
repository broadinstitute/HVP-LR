"""Build manifest of foldseek m8 GCS paths for cohort.

Reads a static list of (sample, workflow_id) pairs (the 48 Succeeded workflows
from Terra submission 96f2f736-f477-4a62-baab-d43158b9278b) and verifies that
each `call-t_08_Format/<sample>.vs_bfvd.m8` blob exists. Writes manifest.tsv.

This avoids 48 sequential Terra API calls; the m8 path is deterministic from
the (submission_id, workflow_id, sample) triple — t_08_Format does not retry
in practice for this pipeline.

History:
- 2026-06-22: initial 45 entries (workflows _4P, _17P, _22P still Running).
- 2026-06-25: added _4P, _17P, _22P after submission completed.
"""
from __future__ import annotations
import sys
from pathlib import Path

from google.cloud import storage

BUCKET = "fc-e71e28f1-164e-4657-a1f3-e7370bad0a0b"
SUBMISSION = "96f2f736-f477-4a62-baab-d43158b9278b"
PREFIX = f"submissions/{SUBMISSION}/HvpViralProteinAnnotation"

# Hard-coded from get_submission_status (all 48 Succeeded as of 2026-06-25).
SAMPLES = [
    ("HVP-0006.1_1P",  "9ec0932c-62ea-41c3-ac2f-292b9928e0e4"),
    ("HVP-0006.1_2P",  "e4bd2497-f5b3-40a4-aacb-422489a8bb59"),
    ("HVP-0006.1_3P",  "f50f8614-02fb-409e-b5a2-e825d2a43997"),
    ("HVP-0006.1_4P",  "9c1174fb-78af-4c31-8505-e6034bbb4e05"),
    ("HVP-0006.1_5P",  "0fb94758-632f-4e03-8ed7-f90d591ed396"),
    ("HVP-0006.1_6P",  "ef5c1776-17d6-4b31-9722-cdc127a08a75"),
    ("HVP-0006.1_7P",  "b332fa15-fbf7-4e84-a92a-9cb86839f191"),
    ("HVP-0006.1_8P",  "0aebd8fc-e0ef-446a-8484-fc74b9a8bed5"),
    ("HVP-0006.1_9P",  "189c00f5-7744-4e16-a43b-cda539a51fb4"),
    ("HVP-0006.1_10P", "78ac5857-5ef6-4fde-99e7-3ad1b9e60dd4"),
    ("HVP-0006.1_11P", "75ddd555-fda7-4281-a35d-3aedc6fed140"),
    ("HVP-0006.1_12P", "99043861-4208-4289-9b8d-9d4bad621ce1"),
    ("HVP-0006.1_13P", "6c21ac3a-73ba-43b3-b6a1-ac98a63f28bb"),
    ("HVP-0006.1_14P", "9876040b-db80-4c7b-b248-4da07e61c555"),
    ("HVP-0006.1_15P", "736ce7a5-d526-4ada-a57d-ac67a5a9a3a7"),
    ("HVP-0006.1_16P", "bec86907-20ef-4871-b157-b24bdd47f65d"),
    ("HVP-0006.1_17P", "e1332a54-8e9e-4eb4-892a-548d6f61abc9"),
    ("HVP-0006.1_18P", "20373acb-c607-4356-a5bc-3d0eed6ea251"),
    ("HVP-0006.1_19P", "2694c304-9729-4313-b67b-ab9b64423457"),
    ("HVP-0006.1_20P", "c63e177d-2776-40a1-878c-e055d29214b7"),
    ("HVP-0006.1_21P", "24879ee9-1a7f-4d9b-9f48-f8383bc0b113"),
    ("HVP-0006.1_22P", "cfefddad-8fb3-4f7b-be00-c6548b15dfc1"),
    ("HVP-0006.1_23P", "8233cc1e-4b9a-4722-85e2-df42074b728c"),
    ("HVP-0006.1_24P", "ecb3e932-0b9b-44d8-a1bc-d5ed71ecc86c"),
    ("HVP-0006.1_25P", "56736622-9ae0-4ece-90f6-229b5d9b55f4"),
    ("HVP-0006.1_26P", "d45fced1-514f-4c29-a0bb-3c01880897c6"),
    ("HVP-0006.1_27P", "fe586d4d-dfec-4f9d-90c0-dbd72680d4ac"),
    ("HVP-0006.1_28P", "fe9d6dd8-3128-40c1-8dfa-c4a65c46abed"),
    ("HVP-0006.1_29P", "9e69299e-3643-4ba5-a1f6-e4dbca8c033c"),
    ("HVP-0006.1_30P", "c40014be-1ce3-4fa0-a523-c53a052b3499"),
    ("HVP-0006.1_31P", "8b581b65-ec49-4a86-85c3-9f325af3d83d"),
    ("HVP-0006.1_32P", "4cfb664d-da31-4f53-86d1-641fc222e660"),
    ("HVP-0006.1_33P", "100b5b1d-5a37-49ac-838f-a9fe672f9295"),
    ("HVP-0006.1_34P", "629d635a-6bc6-489d-849e-d49485c59cfb"),
    ("HVP-0006.1_35P", "e0ec5d0d-142a-4831-bae1-4c041bdaa7f2"),
    ("HVP-0006.1_36P", "6e13cfb7-b0d0-4d41-b713-87014821c0ea"),
    ("HVP-0006.1_37P", "2781d341-30ea-4ddd-bc0e-e6e43d77b911"),
    ("HVP-0006.1_38P", "80b4ca64-ae3d-4aa3-b46f-b124c5f2c1a1"),
    ("HVP-0006.1_39P", "b40e8e08-45fe-4eb4-8da9-69786147184f"),
    ("HVP-0006.1_40P", "2c6e2b63-ffe6-4a51-a70f-eeabe23eaf10"),
    ("HVP-0006.1_41P", "c9be4688-baaa-447a-897c-33fd3eae47a2"),
    ("HVP-0006.1_42P", "ec042ffe-cc34-4a0d-96cd-3668167104f8"),
    ("HVP-0006.1_43P", "8ca76695-fa78-491e-a591-fdc315d3e104"),
    ("HVP-0006.1_44P", "88ecd9c1-1460-4ceb-a177-7524c81c7a1b"),
    ("HVP-0006.1_45P", "6e3538f6-91ef-46f9-9fd9-5bb86520574c"),
    ("HVP-0006.1_46P", "c05801c3-6cd8-4164-93b1-79e310305e2f"),
    ("HVP-0006.1_47P", "d004112f-4826-4034-9d64-676baa711be9"),
    ("HVP-0006.1_48P", "ba7f5593-d305-441c-b619-a0cb30acc388"),
]


def main(out_tsv: Path) -> None:
    client = storage.Client(project="broad-dsde-methods")
    bkt = client.bucket(BUCKET)

    rows = []
    missing = []
    for sample, wfid in SAMPLES:
        blob_path = f"{PREFIX}/{wfid}/call-t_08_Format/{sample}.vs_bfvd.m8"
        blob = bkt.blob(blob_path)
        # blob.exists() is one HEAD request — cheap.
        exists = blob.exists()
        size = blob.size if exists else None
        if not exists:
            missing.append((sample, wfid))
            print(f"[manifest] MISSING {sample} {blob_path}", file=sys.stderr)
            continue
        rows.append((sample, wfid, f"gs://{BUCKET}/{blob_path}", size))
        print(f"[manifest] {sample}  {size/1e6:>7.1f} MB", file=sys.stderr)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as fh:
        fh.write("sample\tworkflow_id\tm8_uri\tsize_bytes\n")
        for sample, wfid, uri, size in rows:
            fh.write(f"{sample}\t{wfid}\t{uri}\t{size}\n")
    total = sum(r[3] for r in rows) / 1e9
    print(f"\n[manifest] wrote {len(rows)}/{len(SAMPLES)} rows  total {total:.2f} GB", file=sys.stderr)
    if missing:
        print(f"[manifest] {len(missing)} missing: {missing}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    proj_root = here.parent
    main(proj_root / "data" / "cohort_2026-06-22" / "manifest.tsv")
