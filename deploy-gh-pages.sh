#!/usr/bin/env bash
# Build and publish docs without changing the source checkout or rewriting history.
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
python_bin=${PYTHON:-python}
python_bin=$(command -v "$python_bin")
export PATH="$(dirname "$python_bin"):$PATH"
build_dir=$(mktemp -d)
site_dir=$(mktemp -d)
cleanup() {
    git worktree remove --force "$site_dir/site" >/dev/null 2>&1 || true
    rm -rf "$build_dir" "$site_dir"
}
trap cleanup EXIT

"$python_bin" -m sphinx -W --keep-going -b html docs "$build_dir/html"
test -s "$build_dir/html/index.html"
git fetch origin gh-pages
git worktree add --detach "$site_dir/site" origin/gh-pages
mkdir -p "$site_dir/site/docs"
rm -rf "$site_dir/site/docs/html"
cp -R "$build_dir/html" "$site_dir/site/docs/html"
touch "$site_dir/site/.nojekyll"
cat > "$site_dir/site/index.html" <<'HTML'
<!doctype html>
<html lang="en"><meta charset="utf-8"><title>Varlens documentation</title>
<meta http-equiv="refresh" content="0; url=docs/html/">
<a href="docs/html/">Varlens documentation</a></html>
HTML
git -C "$site_dir/site" add .
if ! git -C "$site_dir/site" diff --cached --quiet; then
    revision=$(git rev-parse --short HEAD)
    git -C "$site_dir/site" commit -m "Publish Varlens documentation from $revision"
    git -C "$site_dir/site" push origin HEAD:gh-pages
fi
