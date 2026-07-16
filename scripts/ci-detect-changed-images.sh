#!/usr/bin/env bash
# Detect which docker/<image>/ directories changed between a base ref and HEAD.
# Emits to $GITHUB_OUTPUT:
#   images=<json-array>         e.g. ["hvp-monolith","foo"]
#   any=true|false
#   build_matrix=<json-array>   per-(image,arch) build records, e.g.
#       [{"image":"foo","arch":"amd64","runner":"ubuntu-latest","platform":"linux/amd64"},
#        {"image":"foo","arch":"arm64","runner":"ubuntu-24.04-arm","platform":"linux/arm64"}]
#     Drives the build-per-arch matrix in docker.yml. Per-image platform list
#     is read from `PLATFORMS = ...` in each image's Makefile; default both
#     arches if absent or unparseable.
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

# Map a `linux/<arch>` platform string to (arch-name, GH-Actions runner).
# Stdout: "<arch> <runner>"; empty if unknown.
platform_runner() {
    case "$1" in
        linux/amd64) echo "amd64 ubuntu-latest" ;;
        linux/arm64) echo "arm64 ubuntu-24.04-arm" ;;
        *) echo "" ;;
    esac
}

# Build the per-(image, arch) matrix JSON from each image's Makefile
# PLATFORMS line. Default both arches when the line is absent / unparseable
# so existing images keep current behavior with no Makefile edit.
build_matrix_json() {
    local _images_ref="$1"
    eval "local arr=(\"\${${_images_ref}[@]}\")"
    local json="["
    local first=1
    local img mk plats p map arch runner
    for img in "${arr[@]}"; do
        [[ -n "$img" ]] || continue
        mk="docker/$img/Makefile"
        if [[ -f "$mk" ]]; then
            plats=$(awk -F'=' '/^[[:space:]]*PLATFORMS[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$mk")
        else
            plats=""
        fi
        [[ -z "$plats" ]] && plats="linux/amd64,linux/arm64"
        IFS=',' read -ra pa <<<"$plats"
        for p in "${pa[@]}"; do
            map="$(platform_runner "$p")"
            if [[ -z "$map" ]]; then
                echo "warning: $img PLATFORMS entry '$p' is not linux/amd64 or linux/arm64; skipping" >&2
                continue
            fi
            arch="${map% *}"
            runner="${map#* }"
            (( first )) || json+=","
            json+="{\"image\":\"$img\",\"arch\":\"$arch\",\"runner\":\"$runner\",\"platform\":\"$p\"}"
            first=0
        done
    done
    json+="]"
    printf '%s' "$json"
}

emit() {
    local any="$1"
    local images_json="$2"
    local matrix_json="$3"
    if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
        {
            echo "any=$any"
            echo "images=$images_json"
            echo "build_matrix=$matrix_json"
        } >> "$GITHUB_OUTPUT"
    fi
    echo "any=$any"
    echo "images=$images_json"
    echo "build_matrix=$matrix_json"
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
    emit "false" "[]" "[]"
    exit 0
fi

# Build JSON array without depending on jq (avoid adding a runtime dep)
json="["
for i in "${!images[@]}"; do
    [[ $i -gt 0 ]] && json+=","
    json+="\"${images[$i]}\""
done
json+="]"

matrix="$(build_matrix_json images)"

emit "true" "$json" "$matrix"
