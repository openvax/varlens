#!/usr/bin/env bash
# Install this checkout in editable mode with the test and fusion extras.
# Override the interpreter with PYTHON=/path/to/python ./develop.sh
set -euo pipefail

cd "$(dirname "$0")"
"${PYTHON:-python}" -m pip install -e '.[test,fusion]'
