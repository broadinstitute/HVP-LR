"""Integration smoke test for the SexChromKaryotype WDL + x_dosage image.

Runs the full workflow (mosdepth -> rx_sex classify) on a tiny synthetic
Illumina-like BAM whose X/Y read dosage encodes a 46,XY sample, and asserts the
pipeline recovers 46,XY with high confidence. This exercises the real WDL wiring
(both docker images, the TSV->typed-output parse) end-to-end, not just the
classifier unit tests in docker/x_dosage/tests/.

Requirements to actually run (otherwise the test SKIPS, so it is CI-safe):
  - docker available,
  - miniwdl available,
  - the x_dosage image present locally at the tag the WDL references
    (build it: `cd docker/x_dosage && make build`, then
     `docker tag x_dosage:<VER> <GAR_REF>`),
  - network to pull the pinned public mosdepth image (or it present locally).

Fixture: tests/wdl/data/synth_xy.bam(.bai) — see data/make_synth_bam.py for how
it was generated (full-span autosome:X:Y depth 4:2:2 on chr1-3/chrX/chrY).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WDL = os.path.join(REPO_ROOT, "wdl/pipelines/TechAgnostic/QC/SexChromKaryotype.wdl")
BAM = os.path.join(os.path.dirname(__file__), "data/synth_xy.bam")
BAI = os.path.join(os.path.dirname(__file__), "data/synth_xy.bam.bai")
X_DOSAGE_IMAGE = ("us-central1-docker.pkg.dev/broad-hvp-dasc/"
                  "hvp-longread-containers/x_dosage:0.1.0")


def _have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _image_present(ref: str) -> bool:
    try:
        r = subprocess.run(["docker", "image", "inspect", ref],
                           capture_output=True, timeout=30)
        return r.returncode == 0
    except Exception:
        return False


requirements = pytest.mark.skipif(
    not (_have("docker") and _have("miniwdl") and os.path.exists(BAM)
         and _image_present(X_DOSAGE_IMAGE)),
    reason="needs docker + miniwdl + the x_dosage image built/tagged locally "
           f"({X_DOSAGE_IMAGE}); see module docstring",
)


@requirements
def test_workflow_calls_46_xy():
    with tempfile.TemporaryDirectory() as d:
        os.chmod(d, 0o777)  # miniwdl's non-root task user must write here
        proc = subprocess.run(
            ["miniwdl", "run", WDL,
             f"aligned_bam={BAM}", f"aligned_bai={BAI}",
             "sample_name=SYNTH_XY", "--dir", d],
            capture_output=True, text=True, timeout=1200,
        )
        assert proc.returncode == 0, f"miniwdl run failed:\n{proc.stderr[-3000:]}"

        outputs = json.loads(proc.stdout)["outputs"]
        call = outputs["SexChromKaryotype.karyotype_call"]
        conf = outputs["SexChromKaryotype.confidence"]
        rx = outputs["SexChromKaryotype.rx"]
        ry = outputs["SexChromKaryotype.ry"]

        assert call == "46,XY", f"got {call} (rx={rx}, ry={ry})"
        assert conf > 0.9, f"low confidence {conf}"
        assert 0.4 < rx < 0.6, f"Rx off expectation: {rx}"
        assert 0.4 < ry < 0.6, f"Ry off expectation: {ry}"
