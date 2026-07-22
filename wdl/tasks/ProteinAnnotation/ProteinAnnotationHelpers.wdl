version 1.0

import "../../structs/Structs.wdl"

task ConcatProteinFastas {

    meta {
        description: "Concatenate multiple amino-acid FASTA files (any may be .gz) into a single FASTA. Each input is associated with a source label (e.g. vs2, assembly, rescued, genomad); the label is prepended to every sequence header so downstream tools can attribute hits back to the source after clustering. Source labels and FASTAs must be arrays of equal length. Empty inputs are silently skipped."

        tool:         "awk"
        tool_version: "gnu"

        outputs: {
            combined_faa: "Single amino-acid FASTA with source-prefixed headers, ready for mmseqs clustering or foldseek DB creation.",
            num_inputs:   "Number of non-empty input FASTAs concatenated.",
            num_seqs:     "Total sequence count in combined_faa."
        }
    }

    parameter_meta {
        protein_fastas:        "Per-source amino-acid FASTAs to concatenate (any may be .gz)."
        source_labels:         "Per-FASTA source label prepended to each header (e.g. 'vs2', 'genomad', 'assembly', 'rescued'). Must be the same length as protein_fastas."
        prefix:                "Basename for the combined output FASTA."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File]   protein_fastas
        Array[String] source_labels
        String        prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(protein_fastas, "GB"))

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

        OUT=~{prefix}.combined.faa
        rm -f "${OUT}"
        touch "${OUT}"

        num_inputs=0

        # Write files to disk to iterate over:
        protein_fasta_file=~{write_lines(protein_fastas)}
        source_labels_file=~{write_lines(source_labels)}

        # Loop by index using simple variable expansion
        i=0
        while [[ $i -lt ~{length(protein_fastas)} ]]; do
            # Inputs mapping to shell positional parameter access
            FA=$(sed -n "$((i+1))p" "${protein_fasta_file}")
            LB=$(sed -n "$((i+1))p" "${source_labels_file}")

            # Remove trailing .gz if present.
            FA_BASENAME="${FA}"
            if [[ "${FA}" == *.gz ]]; then
                FA_BASENAME="${FA%.gz}"
                gunzip "${FA}"
            fi

            # Now check for records in the uncompressed file, then do annotation
            HAS_RECS=$(grep -m 1 '^>' "$FA_BASENAME")
            if [[ -n "$HAS_RECS" ]]; then
                awk -v L="$LB" '/^>/ { sub(/^>/, ">" L "|", $0) } { print }' "$FA_BASENAME" >> "$OUT"
                num_inputs=$((num_inputs+1))
            else
                echo "INFO: input $FA contains no FASTA records, skipping" >&2
            fi

            # Remove temp file if we created it
            if [[ "$FA_BASENAME" != "$FA" ]]; then
                rm -f "$FA_BASENAME"
            fi
       
            i=$((i+1))
        done

        echo "${num_inputs}" > num_inputs.txt
        awk '/^>/{n++} END{print n+0}' "${OUT}" > num_seqs.txt
   
    >>>

    output {
        File combined_faa = "~{prefix}.combined.faa"
        Int  num_inputs   = read_int("num_inputs.txt")
        Int  num_seqs     = read_int("num_seqs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
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


task AnnotationTransfer {

    meta {
        description: "Reduce a foldseek alignment-hits TSV to one best hit per query (lowest e-value, ties broken by higher bit score) and, when a reference metadata TSV is supplied, left-join the best hit's target ID against that metadata to attach descriptive annotations (organism, function, etc.). The foldseek hits TSV must be headerless and tab-delimited; column names are passed via the hits_columns input and must match the FoldseekConvertAlis format_output used to produce the file (the first two columns must be query and target, in that order, and there must be at least an evalue column for the best-hit reduction; bits is used as the tie-breaker when present). The reference TSV must be tab-delimited with a header row whose first column is the target ID used as the join key. When no reference TSV is given, the annotated output is identical to best_hits_tsv."

        tool:         "python3"
        tool_version: "3.x"

        outputs: {
            best_hits_tsv:       "One row per query — best foldseek hit by lowest e-value. Has a header row matching hits_columns.",
            annotated_hits_tsv:  "best_hits_tsv left-joined against the reference metadata TSV on the target ID. If no reference TSV was provided, identical schema to best_hits_tsv.",
            num_queries_hit:     "Distinct query sequences with at least one hit (= number of rows in best_hits_tsv).",
            num_queries_annotated: "Number of best_hits rows whose target matched a row in the reference metadata TSV. Zero if no reference TSV was provided."
        }
    }

    parameter_meta {
        foldseek_hits_tsv:              "Foldseek alignment hits TSV (output of FoldseekConvertAlis). Headerless; column order must match hits_columns."
        hits_columns:                   "Column names matching the FoldseekConvertAlis format_output of foldseek_hits_tsv. First two entries must be 'query' and 'target'. Must include 'evalue'; if 'bits' is present, it is used as the tie-breaker. Default: query,target,evalue,bits,fident,alnlen,prob (the FoldseekConvertAlis default)."
        reference_metadata_tsv:         "Optional tab-delimited reference metadata TSV; first column is the target ID. Used to left-join annotations onto best hits. Leave unset to skip the join step."
        reference_metadata_has_header:  "Whether reference_metadata_tsv has a header row. When true (default), the first row supplies the metadata column names. When false, supply reference_metadata_columns explicitly (or accept auto-generated meta_col_<N> names). The first row is treated as data when false."
        reference_metadata_columns:     "Optional explicit column names for reference_metadata_tsv. Used when reference_metadata_has_header is false. If empty and no header is provided, synthetic names ('meta_col_1', 'meta_col_2', ...) are generated from the first data row's column count."
        prefix:                         "Basename for the output TSVs."
        runtime_attr_override:          "Override the default runtime attributes."
    }

    input {
        File          foldseek_hits_tsv
        Array[String] hits_columns                  = ["query", "target", "evalue", "bits", "fident", "alnlen", "prob"]
        File?         reference_metadata_tsv
        Boolean       reference_metadata_has_header = true
        Array[String] reference_metadata_columns    = []
        String        prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * (size(foldseek_hits_tsv, "GB") + size(reference_metadata_tsv, "GB")))

    # Python logic in the command block is mirrored in
    # tests/python/annotation_transfer_core.py. Edits here MUST be
    # reflected there (and vice versa) — see tests/python/test_annotation_transfer.py.
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

        python3 - <<'PYEOF'
import csv, sys, os

HITS_PATH       = "~{foldseek_hits_tsv}"
META_PATH       = "~{default='' reference_metadata_tsv}"
BEST_OUT        = "~{prefix}.best_hits.tsv"
ANN_OUT         = "~{prefix}.annotated_hits.tsv"

HITS_COLS = ["~{sep="\", \"" hits_columns}"]
NCOL = len(HITS_COLS)

# Validate schema: first two columns must be query, target; must contain evalue.
if NCOL < 3 or HITS_COLS[0] != "query" or HITS_COLS[1] != "target":
    sys.exit(f"ERROR: hits_columns must start with ['query', 'target', ...]; got {HITS_COLS}")
if "evalue" not in HITS_COLS:
    sys.exit(f"ERROR: hits_columns must include 'evalue'; got {HITS_COLS}")

EVALUE_IDX = HITS_COLS.index("evalue")
BITS_IDX   = HITS_COLS.index("bits") if "bits" in HITS_COLS else None

# 1. Best hit per query (lowest evalue, ties -> higher bits if column present)
best = {}
with open(HITS_PATH, newline="") as f:
    for row in csv.reader(f, delimiter="\t"):
        if not row or len(row) < NCOL:
            continue
        q = row[0]
        try:
            ev_f = float(row[EVALUE_IDX])
        except ValueError:
            continue
        if BITS_IDX is not None:
            try:
                bt_f = float(row[BITS_IDX])
            except ValueError:
                bt_f = 0.0
            key = (ev_f, -bt_f)
        else:
            key = (ev_f, 0.0)
        if q not in best or key < best[q][0]:
            best[q] = (key, row[:NCOL])

with open(BEST_OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(HITS_COLS)
    for q in sorted(best):
        w.writerow(best[q][1])

n_queries_hit = len(best)
print(f"Wrote {BEST_OUT}: {n_queries_hit} queries with best hit", file=sys.stderr)

# 2. Optional annotation join
META_HAS_HEADER  = ~{true="True" false="False" reference_metadata_has_header}
EXPLICIT_COLS = [c for c in ["~{sep="\", \"" reference_metadata_columns}"] if c]

n_annotated = 0
if META_PATH and os.path.exists(META_PATH) and os.path.getsize(META_PATH) > 0:
    meta = {}
    meta_header = None
    n_meta_rows = 0
    n_dups = 0
    with open(META_PATH, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        if META_HAS_HEADER:
            meta_header = next(reader, None)
            if not meta_header:
                print("WARNING: empty reference metadata, falling back to no-join output", file=sys.stderr)
        for r in reader:
            if not r:
                continue
            n_meta_rows += 1
            if meta_header is None and not EXPLICIT_COLS:
                # First data row -> seed synthetic header
                meta_header = [f"meta_col_{i+1}" for i in range(len(r))]
            if r[0] in meta:
                n_dups += 1
            meta[r[0]] = r
    if not META_HAS_HEADER and EXPLICIT_COLS:
        meta_header = list(EXPLICIT_COLS)
    print(f"Loaded {n_meta_rows} metadata rows; {len(meta)} unique keys; {n_dups} duplicate keys (last-row-wins).", file=sys.stderr)
    if meta_header is None:
        out_header = HITS_COLS
        rows = [best[q][1] for q in sorted(best)]
    else:
        ann_cols = ["target_meta_" + c for c in meta_header]
        out_header = HITS_COLS + ann_cols
        rows = []
        for q in sorted(best):
            row = list(best[q][1])
            target = row[1]
            m = meta.get(target)
            if m is not None:
                # Pad/truncate to ann_cols length
                m_padded = (m + [""] * len(ann_cols))[:len(ann_cols)]
                row += m_padded
                n_annotated += 1
            else:
                row += [""] * len(ann_cols)
            rows.append(row)
else:
    out_header = HITS_COLS
    rows = [best[q][1] for q in sorted(best)]

with open(ANN_OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(out_header)
    for r in rows:
        w.writerow(r)
print(f"Wrote {ANN_OUT}: {len(rows)} rows; {n_annotated} annotated via reference", file=sys.stderr)

with open("num_queries_hit.txt", "w") as f:
    f.write(str(n_queries_hit) + "\n")
with open("num_queries_annotated.txt", "w") as f:
    f.write(str(n_annotated) + "\n")
PYEOF
    >>>

    output {
        File best_hits_tsv         = "~{prefix}.best_hits.tsv"
        File annotated_hits_tsv    = "~{prefix}.annotated_hits.tsv"
        Int  num_queries_hit       = read_int("num_queries_hit.txt")
        Int  num_queries_annotated = read_int("num_queries_annotated.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers/hvp-monolith:0.0.3"
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
