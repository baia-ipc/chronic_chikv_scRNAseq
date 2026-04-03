#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"
cd "${RESULTS_DIR}"

mv 037 037-A
mv S3-A3 037-6m
mv 217 217-A
mv 229 229-A
mv 372 372-A
mv S10-A12 267-6m
mv S11-B1 246-A
mv S12-B2 246-6m
mv S13-B3 227-A
mv S14-B4 227-6m
mv S15-A7 041-A
mv S16-A8 041-6m
mv S1-A1 262-6m
mv S2-A2 262-A
mv S4-A4 229-6m
mv S5-A5 217-6m
mv S6-A6 219-6m
mv S7-A9 266-A
mv S8-A10 266-6m
mv S9-A11 267-A
