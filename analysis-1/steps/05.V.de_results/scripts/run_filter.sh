#!/bin/bash

set -euo pipefail

# batch_filter_deseq.sh
# Run filter_unfiltered_deseq_results.py for matched DESeq2 result files in batch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results/tables"
cd "${RESULTS_DIR}"

# Set paths to the directories containing input files
FILTERED_DIR="p_adj_filt"
UNFILTERED_DIR="unfiltered"
OUTPUT_DIR="filtered_unfiltered"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Path to the Python script
SCRIPT="${SCRIPT_DIR}/filter_unfiltered_deseq_results.py"

# Function to run the Python script with given files
run_filter() {
    local filtered="$1"
    local unfiltered="$2"
    local outname="$3"
    python3 "$SCRIPT" "$FILTERED_DIR/$filtered" "$UNFILTERED_DIR/$unfiltered" \
        --output "$OUTPUT_DIR/$outname"
}

# Get list of all files in filtered dir
all_files=("$FILTERED_DIR"/*.tsv)

for filtered_path in "${all_files[@]}"; do
    fname=$(basename "$filtered_path")
    
    if [[ "$fname" == C.* ]]; then
        match="NC${fname:1}"  # C.* <-> NC.*
        if [[ -f "$UNFILTERED_DIR/$match" ]]; then
            run_filter "$fname" "$match" "${match%.tsv}.significant_in_C.tsv"
        fi
    elif [[ "$fname" == NC.* ]]; then
        match="C${fname:2}"  # NC.* <-> C.*
        if [[ -f "$UNFILTERED_DIR/$match" ]]; then
            run_filter "$fname" "$match" "${match%.tsv}.significant_in_NC.tsv"
        fi
    elif [[ "$fname" == A.* ]]; then
        match="M6${fname:1}"  # A.* <-> M6.*
        if [[ -f "$UNFILTERED_DIR/$match" ]]; then
            run_filter "$fname" "$match" "${match%.tsv}.significant_in_A.tsv"
        fi
    elif [[ "$fname" == M6.* ]]; then
        match="A${fname:2}"  # M6.* <-> A.*
        if [[ -f "$UNFILTERED_DIR/$match" ]]; then
            run_filter "$fname" "$match" "${match%.tsv}.significant_in_M6.tsv"
        fi
    fi
done
