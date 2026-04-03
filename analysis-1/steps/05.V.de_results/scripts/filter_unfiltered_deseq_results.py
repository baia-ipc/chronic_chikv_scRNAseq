#!/usr/bin/env python3
"""
filter_unfiltered_deseq_results.py

Extract lines from an unfiltered DESeq2 result file (File2) that match gene names from a filtered DESeq2 result file (File1), preserving the order of genes in File1.

Usage:
  filter_unfiltered_deseq_results.py [options] <filtered_file> <unfiltered_file>

Arguments:
  <filtered_file>    TSV file with p-value filtered DESeq2 results (must contain 'gene' column)
  <unfiltered_file>  TSV file with unfiltered DESeq2 results (must contain 'gene' column)

Options:
  -o FILE, --output=FILE   Output TSV file [default: stdout]
  -v, --verbose            Enable verbose logging
  -q, --quiet              Suppress non-error messages
  -h, --help               Show this help message and exit
  -V, --version            Show version and exit
"""

import sys
import csv
from docopt import docopt
from loguru import logger

def set_logger(verbose: bool, quiet: bool):
    logger.remove()
    if quiet:
        logger.add(lambda msg: None, level="DEBUG")
    elif verbose:
        logger.add(lambda msg: sys.stderr.write(msg), level="INFO")
    else:
        logger.add(lambda msg: sys.stderr.write(msg), level="WARNING")

def extract_matching_genes_ordered(filtered_path, unfiltered_path, output_path):
    logger.info(f"Reading filtered file: {filtered_path}")
    with open(filtered_path, 'r', newline='') as f1:
        reader = csv.DictReader(f1, delimiter='\t')
        filtered_genes = [row['gene'] for row in reader]
    logger.info(f"Collected {len(filtered_genes)} gene names from filtered file")

    logger.info(f"Reading unfiltered file: {unfiltered_path}")
    with open(unfiltered_path, 'r', newline='') as f2:
        reader = csv.DictReader(f2, delimiter='\t')
        header = reader.fieldnames
        unfiltered_data = {row['gene']: row for row in reader}
    logger.info(f"Unfiltered file contains {len(unfiltered_data)} unique genes")

    matching_rows = [unfiltered_data[gene] for gene in filtered_genes if gene in unfiltered_data]
    logger.info(f"Extracted {len(matching_rows)} matching rows in filtered file order")

    if output_path == 'stdout':
        writer = csv.DictWriter(sys.stdout, fieldnames=header, delimiter='\t')
        writer.writeheader()
        writer.writerows(matching_rows)
    else:
        logger.info(f"Writing output to: {output_path}")
        with open(output_path, 'w', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=header, delimiter='\t')
            writer.writeheader()
            writer.writerows(matching_rows)

if __name__ == '__main__':
    args = docopt(__doc__, version='filter_unfiltered_deseq_results 1.1')
    set_logger(args['--verbose'], args['--quiet'])

    filtered_file = args['<filtered_file>']
    unfiltered_file = args['<unfiltered_file>']
    output_file = args['--output']

    extract_matching_genes_ordered(filtered_file, unfiltered_file, output_file)
