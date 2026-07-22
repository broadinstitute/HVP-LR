version 1.0

import "../../structs/Structs.wdl"

task VirSorter2 {

    meta {
        description: "Run VirSorter2 on a merged assembly+rescue FASTA to identify viral sequences. Handles three known packaging issues in the monolith image: (1) the conda/mamba shim required by Snakemake's conda availability check, (2) the missing 'screed' package in the viral conda env, and (3) sklearn 1.3 incompatibility with VS2 DB models pickled under sklearn 0.22.1 (patched in place before running). The DB is copied to a writable location and deleted after the run to reduce output size. Only dsDNAphage and ssDNA groups are included (RNA/NCLDV/lavidaviridae excluded for runtime)."

        tool:         "VirSorter2"
        tool_version: "2.2.4"
        tool_url:     "https://github.com/jiarong/VirSorter2"
        tool_citation: "Guo J, et al. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome. 2021;9(1):37."

        outputs: {
            score_tsv:       "VirSorter2 final-viral-score.tsv: per-sequence viral scores, group assignments, hallmark counts",
            viral_combined_fa: "VirSorter2 final-viral-combined.fa: viral sequences passing the score threshold; input to CheckV"
        }
    }

    parameter_meta {
        merged_fa_gz:          "Merged assembly+rescue FASTA (gzipped) from HvpReadRescue"
        vs2_db_tgz:            "VirSorter2 database as a single compressed archive (.tar.gz or .tar.zst); extracted at runtime"
        sample_name:           "Sample identifier used as output file prefix"
        min_length:            "Minimum sequence length for VirSorter2 to consider a sequence (default 1000)"
        extra_args:            "Additional command-line args appended verbatim to the virsorter run invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File   merged_fa_gz
        File   vs2_db_tgz
        String sample_name

        Int    min_length = 1000
        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(5.0 * size(vs2_db_tgz, "GB")) + ceil(10.0 * size(merged_fa_gz, "GB"))

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

        # Conda/mamba shim: Snakemake checks for a 'conda' executable; monolith uses micromamba
        mkdir -p /tmp/condabin
        cat > /tmp/condabin/mamba << 'SHIM_EOF'
#!/bin/bash
if [[ "${1:-}" == "--version" ]]; then
    echo "mamba 23.3.1"
else
    exec micromamba "$@"
fi
SHIM_EOF
        cp /tmp/condabin/mamba /tmp/condabin/conda
        chmod +x /tmp/condabin/conda /tmp/condabin/mamba
        export PATH="/tmp/condabin:/opt/conda/envs/viral/bin:$PATH"

        # Fix missing screed and sklearn 1.3 incompatibility with VS2 0.22.1-pickled models
        micromamba install -p /opt/conda/envs/viral -c conda-forge -c bioconda \
            screed "scikit-learn=1.2.*" "imbalanced-learn=0.10.*" --yes --quiet

        # Extract VS2 database and copy to a writable location (required by VS2)
        mkdir -p vs2_db_src
        if [[ "~{vs2_db_tgz}" == *.tar.zst ]]; then
            zstd -d ~{vs2_db_tgz} --stdout | tar -x -C vs2_db_src
        else
            tar -xzf ~{vs2_db_tgz} -C vs2_db_src
        fi
        DB_SRC=$(find vs2_db_src -name 'group' -type d | head -1 | xargs dirname)
        mkdir -p vs2_out
        cp -r "${DB_SRC}" vs2_out/vs2_db
        rm -rf vs2_out/vs2_db/conda_envs

        # Patch VS2 DB models: sklearn 1.3 added missing_go_to_left to DecisionTree node
        # dtype, making 0.22.1-pickled models unloadable. Also patch MinMaxScaler.clip.
        python3 - <<'PYEOF'
import joblib, glob, sys, warnings
from sklearn.preprocessing import MinMaxScaler

def patch_mms(obj):
    if isinstance(obj, MinMaxScaler) and not hasattr(obj, 'clip'):
        obj.clip = False
    if isinstance(obj, dict):
        for v in obj.values():
            patch_mms(v)
    if hasattr(obj, 'steps'):
        for _, step in obj.steps:
            patch_mms(step)
    if hasattr(obj, 'estimators_'):
        for est in obj.estimators_:
            patch_mms(est)

for model_path in glob.glob('vs2_out/vs2_db/group/*/model'):
    print(f'Patching {model_path}...', file=sys.stderr, flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = joblib.load(model_path)
    patch_mms(model)
    joblib.dump(model, model_path)
PYEOF

        # Decompress input FASTA if gzipped — VS2 does grep -c '^>' internally
        CONTIGS_INPUT="~{merged_fa_gz}"
        if [[ "~{merged_fa_gz}" == *.gz ]]; then
            CONTIGS_INPUT="/tmp/input_contigs.fasta"
            gunzip -c ~{merged_fa_gz} > "${CONTIGS_INPUT}"
        fi

        # Dump step3-classify logs to stderr on failure for easier debugging
        trap '
          CLASSIFY_LOGDIR="vs2_out/log/iter-0/step3-classify"
          if [ -d "${CLASSIFY_LOGDIR}" ]; then
            for log in "${CLASSIFY_LOGDIR}"/*.log; do
              [ -f "${log}" ] || continue
              echo "=== $(basename "${log}") ===" >&2
              cat "${log}" >&2
            done
          fi
        ' ERR

        virsorter run \
            -i "${CONTIGS_INPUT}" \
            -w vs2_out \
            -d vs2_out/vs2_db \
            --include-groups dsDNAphage,ssDNA \
            --min-length ~{min_length} \
            --jobs "${NUM_CPUS}" \
            --provirus-off \
            --use-conda-off \
            ~{extra_args}

        rm -rf vs2_out/vs2_db

        cp vs2_out/final-viral-score.tsv    ~{sample_name}.vs2_score.tsv
        cp vs2_out/final-viral-combined.fa  ~{sample_name}.vs2_viral_combined.fa
    >>>

    output {
        File score_tsv        = "~{sample_name}.vs2_score.tsv"
        File viral_combined_fa = "~{sample_name}.vs2_viral_combined.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
