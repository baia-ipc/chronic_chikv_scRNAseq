#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"
cd "${RESULTS_DIR}"

for dir in prev_analysis_1/*/outs; do ln -s "$dir" "$(basename "$(dirname "$dir")")"; done
for dir in prev_analysis_2/*/outs; do ln -s "$dir" "$(basename "$(dirname "$dir")")"; done
