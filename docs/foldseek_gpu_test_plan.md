# Foldseek GPU image — validation plan for a GPU-equipped host

You are a Claude Code instance running on a machine with at least one NVIDIA GPU of compute capability **>= 7.5** (Turing or newer: T4, A100, RTX 30xx/40xx, L4, H100). The CPU-side foldseek image and its WDL surface have already been validated end-to-end on a CPU-only host. This document covers the GPU-only paths that the CPU sandbox could not exercise.

## Prerequisites — verify before testing

1. **NVIDIA driver + CUDA runtime present on host.**

   ```bash
   nvidia-smi
   ```

   Confirm at least one GPU is listed, driver version is >= 545 (CUDA 12.6 minimum), and the GPU is idle. Note the device's compute capability — if it is older than 7.5 the binary will not load and the rest of this plan is moot; switch to a different host.

2. **nvidia-container-toolkit installed.**

   ```bash
   docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu24.04 nvidia-smi
   ```

   Should print the same `nvidia-smi` output from inside the container. If this fails, install `nvidia-container-toolkit` (`apt install nvidia-container-toolkit && nvidia-ctk runtime configure --runtime=docker && systemctl restart docker`) before proceeding.

3. **Image pulled (or built locally).**

   Either pull the published image (after this branch merges):

   ```bash
   docker pull ghcr.io/broadinstitute/hvp-lr/foldseek-gpu:10.0.1
   docker tag ghcr.io/broadinstitute/hvp-lr/foldseek-gpu:10.0.1 foldseek-gpu:10.0.1
   ```

   Or build from the branch:

   ```bash
   cd HVP-LR/docker/foldseek-gpu
   make build
   ```

   Either way, after this step `docker images foldseek-gpu:10.0.1` must show the image at roughly 1.5 GB.

4. **Test fixtures.** A handful of PDB files. The CPU sandbox used:

   ```
   1MBO 2DN1 1HHO 1LH1 2HBG       # globins
   1AKE 1CRN 1UBQ 1LYZ 7TIM       # unrelated folds
   ```

   Fetch them once:

   ```bash
   mkdir -p /tmp/fseek_gpu/pdbs && cd /tmp/fseek_gpu/pdbs
   for id in 1MBO 2DN1 1HHO 1LH1 2HBG 1AKE 1CRN 1UBQ 1LYZ 7TIM; do
       curl -sSL "https://files.rcsb.org/download/${id}.pdb" -o "${id}.pdb"
   done
   ls -l
   ```

## Test 1 — image loads CUDA at runtime

Goal: confirm the binary actually initialises CUDA when --gpus is passed and gracefully refuses (or falls back) when it is not.

```bash
# With GPU:
docker run --rm --gpus all foldseek-gpu:10.0.1 version
# Expect: 718d42176d2f67d36a60866fedfb881f8d5a7ebf

# Without --gpus (CPU-only command):
docker run --rm foldseek-gpu:10.0.1 easy-search /dev/null /dev/null /tmp/x /tmp/y 2>&1 | head -5
# Expect: foldseek prints help / errors on the missing files, NOT a CUDA init error.
# The image must remain usable for CPU-only subcommands even when no GPU is attached.
```

**Pass criteria.** `version` returns the SHA verbatim. CPU subcommands run without `--gpus` and do not crash on missing CUDA. Record the exact output for both invocations.

## Test 2 — ProstT5 download + GPU inference end-to-end

The point of this image. Downloads the ProstT5 weights (~1 GB), then uses GPU-accelerated sequence-to-3Di inference to build a structure DB from an amino-acid FASTA, and searches a real PDB against it.

```bash
mkdir -p /tmp/fseek_gpu/prostt5 /tmp/fseek_gpu/work
cd /tmp/fseek_gpu

# Step A: download ProstT5 model weights via foldseek's databases command.
docker run --rm --gpus all \
    -v /tmp/fseek_gpu:/work \
    foldseek-gpu:10.0.1 \
    databases ProstT5 /work/prostt5/weights /work/prostt5_tmp

# Confirm:
ls -lh /tmp/fseek_gpu/prostt5/
# Expect: weights file(s) totalling ~1 GB.

# Step B: build a query DB from amino-acid sequences using ProstT5 on GPU.
# Extract sequences from the globin PDBs into a FASTA so there is a real
# AA-only query for the ProstT5 path.
docker run --rm \
    -v /tmp/fseek_gpu:/work \
    foldseek-gpu:10.0.1 \
    convert2pdb /work/work/dummy /work/work/scratch 2>&1 || true
# (Easier alternative: ship a small AA fasta directly.)
cat > /tmp/fseek_gpu/work/queries.fasta <<'EOF'
>1MBO_sequence
GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG
>1HHO_alpha
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR
EOF

docker run --rm --gpus all \
    -v /tmp/fseek_gpu:/work \
    foldseek-gpu:10.0.1 \
    createdb /work/work/queries.fasta /work/work/queryDB \
    --prostt5-model /work/prostt5/weights --gpu 1

# Confirm:
ls /tmp/fseek_gpu/work/queryDB*
# Expect: foldseek DB files including queryDB.dbtype, queryDB_ss (3Di sequences
# predicted by ProstT5), queryDB_h, etc.

# Step C: build a target DB from real PDBs (no ProstT5 needed for targets).
docker run --rm --gpus all \
    -v /tmp/fseek_gpu:/work \
    -v /tmp/fseek_gpu/pdbs:/pdbs \
    foldseek-gpu:10.0.1 \
    createdb /pdbs /work/work/targetDB

# Step D: GPU-accelerated structural search.
docker run --rm --gpus all \
    -v /tmp/fseek_gpu:/work \
    foldseek-gpu:10.0.1 \
    easy-search \
        /work/work/queryDB \
        /work/work/targetDB \
        /work/work/result.m8 \
        /work/work/search_tmp \
        --gpu 1 -e 1e-3
```

**Pass criteria.**
- ProstT5 weights download completes and produces a ~1 GB file.
- `createdb --prostt5-model --gpu 1` finishes, GPU is visibly utilised during the call (`nvidia-smi -l 1` in another terminal should show foldseek_avx2 holding the GPU), and the resulting queryDB contains 3Di sequences (check `queryDB_ss` exists and is non-empty).
- `easy-search --gpu 1` returns hits including 1MBO_sequence -> 1MBO and 1HHO_alpha -> 1HHO_A as top hits. E-values must be small (<1e-15 for self-matches if PDB chain matches the source sequence).

**Time budget.** ProstT5 download: ~5 min. ProstT5 inference on two sequences: < 30 s on any supported GPU. Search: a few seconds.

## Test 3 — verify --gpu 0 falls back to CPU code path

Same image, GPU flag off, no `--gpus` runtime flag. Confirms users on the GPU image can still execute CPU-only workflows when they choose to.

```bash
docker run --rm \
    -v /tmp/fseek_gpu/pdbs:/pdbs \
    foldseek-gpu:10.0.1 \
    easy-search /pdbs/1MBO.pdb /pdbs/1HHO.pdb /tmp/result.m8 /tmp/tmp_cpu
```

**Pass criteria.** Same hits as the CPU image would produce (1MBO -> 1HHO_A at evalue ~ 4e-7). Diff against the CPU-image's output of the same command if you have it handy.

## Test 4 — WDL task using the GPU image

Exercise the runtime-override mechanism (option A from the design): a WDL task default-pointing at the CPU image gets redirected to the GPU image via inputs JSON.

```bash
cat > /tmp/wdl_gpu_in.json <<'EOF'
{
  "FoldseekStructuralSearch.query_structures": ["/tmp/fseek_gpu/pdbs/1MBO.pdb"],
  "FoldseekStructuralSearch.target_structures": [
    "/tmp/fseek_gpu/pdbs/1MBO.pdb",
    "/tmp/fseek_gpu/pdbs/2DN1.pdb",
    "/tmp/fseek_gpu/pdbs/1HHO.pdb"
  ],
  "FoldseekStructuralSearch.prefix": "gpu_wdl_smoke",
  "FoldseekStructuralSearch.t_01_FoldseekEasySearch.runtime_attr_override": {
    "docker": "foldseek-gpu:10.0.1"
  }
}
EOF

cd HVP-LR
miniwdl run wdl/pipelines/TechAgnostic/ProteinAnnotation/FoldseekStructuralSearch.wdl \
    -i /tmp/wdl_gpu_in.json -d /tmp/wdl_gpu_run/
```

**Pass criteria.** Workflow completes. `num_hits` matches the same query/target combo run against the CPU image (should be 5 hits: 1MBO self + 4 cross-globin alignments). This confirms the runtime override mechanism is the correct WDL surface for GPU opt-in.

Note: miniwdl does NOT forward `--gpus` automatically. To actually use the GPU under miniwdl you need either a custom docker runtime or a Cromwell config. The miniwdl smoke run still proves the image is callable from WDL; the GPU activation step belongs to whichever orchestrator drives production runs.

## Test 5 — Terra / Cromwell integration (optional, only if you have a GPU runner)

```
runtime_attr_override = {
  "docker": "ghcr.io/broadinstitute/hvp-lr/foldseek-gpu:10.0.1",
  "gpuCount": 1,
  "gpuType": "nvidia-tesla-t4",
  "cpu_cores": 8,
  "mem_gb": 32
}
```

Submit a single FoldseekEasySearch run against a real query/target pair. Check Cromwell logs for GPU allocation and inspect `/proc/driver/nvidia/version` in the task workspace bucket if available.

## Reporting back

When this plan is complete, append a section to `/workspace/JOURNAL.md` summarising:

- GPU model + driver version + CUDA runtime version on the host.
- Exact image SHA tested (from `docker inspect foldseek-gpu:10.0.1 --format='{{.Id}}'`).
- Pass / fail for each test, with the actual output for tests 1, 2, and 4.
- ProstT5 download size and inference wall time.
- Any deviation from the expected hit count or output schema.

If any test fails, capture the full task `stderr.txt` and reference the call directory path; do not delete the workspace until the failure is triaged. Open a follow-up commit on the branch with a fix proposal — do not push without operator approval (matches the standing security posture).
