#!/bin/bash
set -e

# Paths

REPO_ROOT="$(dirname "$(dirname "$(realpath "$0")")")"
NOTEBOOKS_SRC="$REPO_ROOT/examples/notebooks/ETC_notebooks"
NOTEBOOKS_DST="$REPO_ROOT/docs/source/examples/notebooks"

echo "Starting pre-build processes for Sphinx documentation"
echo ""
echo "Repo root resolved to: $REPO_ROOT"
echo "Notebook source: $NOTEBOOKS_SRC"
echo "Notebook destination: $NOTEBOOKS_DST"
echo ""

echo "Updating git submodules..."
if git submodule update --init --recursive; then
    echo "Git submodules updated successfully."
else
    echo "Failed to update git submodules." >&2
    exit 1
fi
echo ""

echo "Remove old notebooks folder if exists..."
if rm -rf "$NOTEBOOKS_DST"; then
    echo "Success"
else
    echo "Failed to remove old notebooks folder." >&2
    exit 1
fi
echo ""

echo "Copying notebooks from submodule..."

find "$NOTEBOOKS_SRC" -name '*.ipynb' | while read -r file; do
    rel_path="${file#$NOTEBOOKS_SRC/}"
    dest_path="$NOTEBOOKS_DST/$rel_path"
    dest_dir="$(dirname "$dest_path")"

    mkdir -p "$dest_dir"

    cp "$file" "$dest_path"
done

echo ""
