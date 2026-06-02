version 1.0

import "../../structs/Structs.wdl"

# -----------------------------------------------------------------------------
# Ported from the Broad Institute's Long Read Pipelines (LRP) repository.
#
#   LRP repo:       https://github.com/broadinstitute/long-read-pipelines
#   WDL source:     wdl/tasks/Utility/PBUtils.wdl  (task SummarizePBI)
#                   https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/Utility/PBUtils.wdl
#   Python source:  docker/lr-pb/compute_pbi_stats.py
#                   https://github.com/broadinstitute/long-read-pipelines/blob/main/docker/lr-pb/compute_pbi_stats.py
#   Original author: Kiran V Garimella (Broad Institute, LRMA group)
#
# The LRP implementation depends on a custom container image
# (us.gcr.io/broad-dsp-lrma/lr-pb) that bakes compute_pbi_stats.py in.
# This port inlines the script verbatim and runs it under python:3.11-slim
# with construct + numpy pip-installed at task start, so HVP-LR has no
# dependency on the lr-pb image. Output names are also adapted to HVP-LR's
# scalar-output convention (one Int/Float per metric) instead of the
# Map[String, Float] surfaced by LRP's SummarizePBI.
# -----------------------------------------------------------------------------

task SummarizeHifiPbi {

    meta {
        description: "Parse a PacBio .pbi index file and emit per-sample polymerase and subread read-length and quality summary metrics. Ported from broadinstitute/long-read-pipelines (wdl/tasks/Utility/PBUtils.wdl :: SummarizePBI; docker/lr-pb/compute_pbi_stats.py). Original author: Kiran V Garimella."

        tool:         "compute_pbi_stats.py (numpy + construct)"
        tool_url:     "https://github.com/broadinstitute/long-read-pipelines/blob/main/docker/lr-pb/compute_pbi_stats.py"

        outputs: {
            pbi_stats_map:     "Two-column TSV (key<TAB>value) of all PBI-derived metrics",
            pbi_reads:         "Number of reads at or above the quality threshold",
            pbi_bases:         "Total bases across all reads",
            pbi_mean_qual:     "Mean per-read Phred quality score (from the PBI readQual field)",
            pbi_median_qual:   "Median per-read Phred quality score (from the PBI readQual field)",
            polymerase_mean:   "Mean polymerase read length (sum of subread lengths within each ZMW)",
            polymerase_median: "Median polymerase read length",
            polymerase_stdev:  "Standard deviation of polymerase read lengths",
            polymerase_n50:    "N50 polymerase read length",
            subread_mean:      "Mean subread length (one entry per PBI record)",
            subread_median:    "Median subread length",
            subread_stdev:     "Standard deviation of subread lengths",
            subread_n50:       "N50 subread length"
        }
    }

    parameter_meta {
        pbi:                   "PacBio .pbi index file for an unaligned HiFi BAM"
        qual_threshold:        "Phred-scale per-read quality threshold; only reads at or above this threshold are included (default 0 = all reads)"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File pbi
        Int  qual_threshold = 0

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 2 * ceil(size(pbi, "GB"))

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        # Install runtime dependencies. construct is pinned to the same major
        # version family that the LRP lr-pb container uses (2.9.x) to keep the
        # parser API in lockstep with the original implementation.
        pip install --no-cache-dir --quiet 'numpy<2' 'construct<2.10'

        # Inline copy of compute_pbi_stats.py from broadinstitute/long-read-pipelines.
        # Source: https://github.com/broadinstitute/long-read-pipelines/blob/main/docker/lr-pb/compute_pbi_stats.py
        # Original author: Kiran V Garimella.
        cat > compute_pbi_stats.py <<'PYEOF'
from __future__ import print_function
import sys
import argparse
import gzip
import math
import numpy
from construct import *


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def to_phred_score(p):
    return 40 if p >= 1.0 else int(-10.0 * math.log10(1.0 - p))


def n50(lengths):
    all_len = sorted(lengths, reverse=True)
    csum = numpy.cumsum(all_len)
    n2 = int(sum(lengths)/2)
    csumn2 = min(csum[csum >= n2])
    ind = numpy.where(csum == csumn2)

    return all_len[int(ind[0])]


def load_index(pbi_file, qual_threshold):
    """
    Load .pbi index data
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need for sharding.
    # More on index format at https://pacbiofileformats.readthedocs.io/en/9.0/PacBioBamIndex.html .

    fmt = Struct(
        # Header
        "magic" / Const(b"PBI\x01"),
        "version_patch" / Int8ul,
        "version_minor" / Int8ul,
        "version_major" / Int8ul,
        "version_empty" / Int8ul,
        "pbi_flags" / Int16ul,
        "n_reads" / Int32ul,
        "reserved" / Padding(18),

        # Basic information section (columnar format)
        "rgId" / Padding(this.n_reads * 4),
        "qStart" / Array(this.n_reads, Int32sl),
        "qEnd" / Array(this.n_reads, Int32sl),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Array(this.n_reads, Float32l),
    )

    n_reads = 0
    polymerase_read_lengths = {}
    subread_lengths = []
    quals = []
    total_bases = 0
    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            if to_phred_score(idx_contents.readQual[j]) >= qual_threshold:
                length = idx_contents.qEnd[j] - idx_contents.qStart[j]

                # Save the polymerase and subread lengths
                if idx_contents.holeNumber[j] not in polymerase_read_lengths:
                    polymerase_read_lengths[idx_contents.holeNumber[j]] = 0

                n_reads += 1
                polymerase_read_lengths[idx_contents.holeNumber[j]] += length
                subread_lengths.append(length)
                quals.append(to_phred_score(idx_contents.readQual[j]))
                total_bases += length

    return n_reads, total_bases, numpy.mean(quals), numpy.median(quals), polymerase_read_lengths, subread_lengths


def main():
    parser = argparse.ArgumentParser(description='Compute .pbi stats', prog='compute_pbi_stats')
    parser.add_argument('-q', '--qual-threshold', type=int, default=0, help="Phred-scale quality threshold")
    parser.add_argument('pbi', type=str, help=".pbi index")
    args = parser.parse_args()

    # Decode PacBio .pbi file and determine the polymerase and subread lengths
    eprint(f"Reading index ({args.pbi}). This may take a few minutes...", flush=True)
    n_reads, n_bases, mean_qual, median_qual, polymerase_read_lengths, subread_lengths = load_index(args.pbi, args.qual_threshold)

    prl = list(polymerase_read_lengths.values())

    print(f'reads\t{n_reads}')
    print(f'bases\t{n_bases}')
    print(f'mean_qual\t{mean_qual if len(prl) else 0}')
    print(f'median_qual\t{median_qual if len(prl) else 0}')

    print(f'polymerase_mean\t{int(numpy.mean(prl)) if len(prl) else 0}')
    print(f'polymerase_median\t{int(numpy.median(prl)) if len(prl) else 0}')
    print(f'polymerase_stdev\t{int(numpy.std(prl)) if len(prl) else 0}')
    print(f'polymerase_n50\t{n50(prl) if len(prl) else 0}')

    print(f'subread_mean\t{int(numpy.mean(subread_lengths)) if len(subread_lengths) else 0}')
    print(f'subread_median\t{int(numpy.median(subread_lengths)) if len(subread_lengths) else 0}')
    print(f'subread_stdev\t{int(numpy.std(subread_lengths)) if len(subread_lengths) else 0}')
    print(f'subread_n50\t{n50(subread_lengths) if len(subread_lengths) else 0}')


if __name__ == "__main__":
    main()
PYEOF

        python3 compute_pbi_stats.py -q ~{qual_threshold} ~{pbi} | tee pbi_stats.map.txt

        # Split the key<TAB>value map into per-key scalar files for WDL bindings.
        while IFS=$'\t' read -r key val; do
            printf '%s\n' "${val}" > "stat.${key}.txt"
        done < pbi_stats.map.txt
    >>>

    output {
        File pbi_stats_map = "pbi_stats.map.txt"

        Int   pbi_reads         = read_int("stat.reads.txt")
        Int   pbi_bases         = read_int("stat.bases.txt")
        Float pbi_mean_qual     = read_float("stat.mean_qual.txt")
        Float pbi_median_qual   = read_float("stat.median_qual.txt")

        Int polymerase_mean   = read_int("stat.polymerase_mean.txt")
        Int polymerase_median = read_int("stat.polymerase_median.txt")
        Int polymerase_stdev  = read_int("stat.polymerase_stdev.txt")
        Int polymerase_n50    = read_int("stat.polymerase_n50.txt")

        Int subread_mean      = read_int("stat.subread_mean.txt")
        Int subread_median    = read_int("stat.subread_median.txt")
        Int subread_stdev     = read_int("stat.subread_stdev.txt")
        Int subread_n50       = read_int("stat.subread_n50.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "python:3.11-slim"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
