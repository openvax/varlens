#!/bin/bash
set -euo pipefail
python -m ruff check varlens tests setup.py
