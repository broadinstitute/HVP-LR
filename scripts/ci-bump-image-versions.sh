#!/usr/bin/env bash
# Patch-bump VERSION in docker/<image>/Makefile for each image in the given JSON
# array, UNLESS the dev already changed that VERSION line in the pushed range.
#
# Args:
#   $1 : JSON array of image names, e.g. ["hvp-monolith","foo"]
#
# Env:
#   BEFORE_SHA : pre-push SHA for diff (defaults to HEAD~1)
#
# Emits to $GITHUB_OUTPUT:
#   committed=true|false        whether any Makefile was modified
#   versions=<json-object>      {"<image>": "<new-or-current-version>"}
#   bumped=<json-array>         images we actually bumped
set -euo pipefail

images_json="${1:?usage: $0 <json-image-array>}"
base="${BEFORE_SHA:-HEAD~1}"
zero_sha="0000000000000000000000000000000000000000"
if [[ "$base" == "$zero_sha" ]]; then
    base="HEAD~1"
fi

# Parse JSON array without jq dependency: strip brackets/quotes/whitespace and split on ,
images_csv="$images_json"
images_csv="${images_csv//\[/}"
images_csv="${images_csv//\]/}"
images_csv="${images_csv//\"/}"
images_csv="${images_csv// /}"
IFS=',' read -r -a images <<<"$images_csv"

versions_json="{"
bumped_json="["
committed=false
first_ver=true
first_bump=true

for img in "${images[@]}"; do
    [[ -n "$img" ]] || continue
    mk="docker/$img/Makefile"
    if [[ ! -f "$mk" ]]; then
        echo "warning: $mk not found; skipping $img" >&2
        continue
    fi

    cur="$(awk -F'=' '/^[[:space:]]*VERSION[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$mk")"
    if [[ -z "$cur" ]]; then
        echo "warning: no VERSION line in $mk; skipping $img" >&2
        continue
    fi

    # Has the VERSION line itself moved in the pushed range? If so, dev (or a
    # prior bump commit) handled it — leave alone.
    if git rev-parse --verify "$base" >/dev/null 2>&1 \
       && git diff "$base" HEAD -- "$mk" \
            | grep -E '^[-+][[:space:]]*VERSION[[:space:]]*=' >/dev/null 2>&1; then
        echo "$img: VERSION line already changed in range; no auto-bump (current=$cur)"
        new="$cur"
    else
        IFS='.' read -r ma mi pa <<<"$cur"
        if [[ -z "${ma:-}" || -z "${mi:-}" || -z "${pa:-}" ]]; then
            echo "warning: $img VERSION '$cur' is not X.Y.Z; skipping bump" >&2
            new="$cur"
        else
            pa=$((pa + 1))
            new="${ma}.${mi}.${pa}"
            # Portable in-place sed (GNU; CI runs on ubuntu-latest)
            sed -i -E "s|^([[:space:]]*VERSION[[:space:]]*=[[:space:]]*).*$|\1${new}|" "$mk"
            echo "$img: bumped $cur -> $new"
            committed=true
            $first_bump || bumped_json+=","
            bumped_json+="\"$img\""
            first_bump=false
        fi
    fi

    $first_ver || versions_json+=","
    versions_json+="\"$img\":\"$new\""
    first_ver=false
done

versions_json+="}"
bumped_json+="]"

if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
    {
        echo "committed=$committed"
        echo "versions=$versions_json"
        echo "bumped=$bumped_json"
    } >> "$GITHUB_OUTPUT"
fi

echo "committed=$committed"
echo "versions=$versions_json"
echo "bumped=$bumped_json"
