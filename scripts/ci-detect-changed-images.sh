#!/usr/bin/env bash
# Detect which docker/<image>/ directories changed between a base ref and HEAD.
# Emits to $GITHUB_OUTPUT:
#   images=<json-array>   e.g. ["hvp-monolith","foo"]
#   any=true|false
#
# Driven by env vars set by the workflow:
#   GITHUB_EVENT_NAME : push | pull_request | schedule | workflow_dispatch
#   BASE_REF          : PR base ref (for pull_request)
#   BEFORE_SHA        : pre-push SHA (for push)
#
# For schedule and workflow_dispatch we return every docker/*/ dir that has a
# Dockerfile (the scan path uses this to rescan all latest images).
set -euo pipefail

event="${GITHUB_EVENT_NAME:-}"
zero_sha="0000000000000000000000000000000000000000"

list_all_image_dirs() {
    local d
    for d in docker/*/; do
        [[ -f "${d}Dockerfile" ]] || continue
        basename "$d"
    done
}

resolve_base() {
    case "$event" in
        pull_request)
            if [[ -n "${BASE_REF:-}" ]]; then
                git fetch --no-tags --depth=200 origin "$BASE_REF" >/dev/null 2>&1 || true
                echo "origin/${BASE_REF}"
            else
                echo "HEAD~1"
            fi
            ;;
        push)
            if [[ -n "${BEFORE_SHA:-}" && "$BEFORE_SHA" != "$zero_sha" ]]; then
                echo "$BEFORE_SHA"
            else
                # New branch push: compare against parent, falling back to HEAD~1
                if git rev-parse HEAD~1 >/dev/null 2>&1; then
                    echo "HEAD~1"
                else
                    echo ""  # Initial commit; treat as "all images"
                fi
            fi
            ;;
        *)
            echo ""
            ;;
    esac
}

emit() {
    local any="$1"
    local images_json="$2"
    if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
        {
            echo "any=$any"
            echo "images=$images_json"
        } >> "$GITHUB_OUTPUT"
    fi
    echo "any=$any"
    echo "images=$images_json"
}

# Schedule / dispatch / initial-commit push: scan every image dir.
if [[ "$event" == "schedule" || "$event" == "workflow_dispatch" ]]; then
    mapfile -t images < <(list_all_image_dirs)
else
    base="$(resolve_base)"
    if [[ -z "$base" ]]; then
        mapfile -t images < <(list_all_image_dirs)
    else
        # Top-level subdir name under docker/ for any changed file
        mapfile -t changed < <(git diff --name-only "$base" HEAD -- 'docker/*/*' \
            | awk -F/ 'NF>=3 {print $2}' | sort -u)
        images=()
        for img in "${changed[@]}"; do
            [[ -n "$img" && -f "docker/$img/Dockerfile" ]] && images+=("$img")
        done
    fi
fi

if (( ${#images[@]} == 0 )); then
    emit "false" "[]"
    exit 0
fi

# Build JSON array without depending on jq (avoid adding a runtime dep)
json="["
for i in "${!images[@]}"; do
    [[ $i -gt 0 ]] && json+=","
    json+="\"${images[$i]}\""
done
json+="]"

emit "true" "$json"
