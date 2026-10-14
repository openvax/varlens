#!/bin/bash
set -euo pipefail
python -m pytest --cov=varlens --cov-report=term-missing "$@"
