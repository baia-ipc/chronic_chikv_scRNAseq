#!/usr/bin/env python3
"""
Run cellranger count for the samples described in the project metadata table.

For each selected sample one cellranger run directory named after the sample's
SampleID is created under <outdir>; its outs/ subdirectory holds the results.
The FASTQ location of a sample is derived from the metadata columns as:

    <readsrootdir>/<RunID>/<RunID>_10X_RawData_Outs/<DemultiplexedID>/<FlowcellID>

This single convention covers both the 2023 and the 2025 sequencing batches.
Samples whose output directory already exists are skipped, so that previously
computed (and expensive) results are never overwritten unless --force is given.

Usage:
    run_cellranger_count.py [options] \
            <outdir> <metadata_table> <transcrdir> <readsrootdir>

Arguments:
    <outdir>            Output root directory (one subdir per sample is created)
    <metadata_table>    TSV file containing the samples metadata
    <transcrdir>        Cellranger reference transcriptome directory
    <readsrootdir>      Reads root directory

Options:
    -s, --samples LIST     Select samples by SampleID in ,-sep LIST
    -p, --patients LIST    Select samples by PatientID in ,-sep LIST
    -t, --timepoints LIST  Select samples by TimePoint in ,-sep LIST
    -y, --years LIST       Select samples by SequencingYear in ,-sep LIST
    -r, --runs LIST        Select samples by RunID in ,-sep LIST
    --localcores N         Cores passed to cellranger [default: 60]
    --localvmem N          Memory (GB) passed to cellranger [default: 200]
    -f, --force            Run even if the output directory already exists
    -D, --dry              Dry run (use --dry option in cellranger)
    -v, --verbose          Print more output
    -d, --debug            Print debug output
    -q, --quiet            Print less output
    -h, --help             Show this screen
    -V, --version          Show version
"""

import csv
import docopt
import sh
import sys
import re
import os
from loguru import logger


def setup_logger(args):
    if sum([args['--debug'], args['--verbose'], args['--quiet']]) > 1:
        logger.error("Only one of --debug, --verbose, --quiet allowed")
        sys.exit(1)
    logger.remove()
    loglevel = 'WARNING'
    if args['--debug']: loglevel = 'DEBUG'
    if args['--verbose']: loglevel = 'INFO'
    if args['--quiet']: loglevel = 'ERROR'
    logger.add(sys.stderr, level=loglevel)


def postprocess_args(args):
    for pathkey in ['<outdir>', '<metadata_table>', '<transcrdir>',
                    '<readsrootdir>']:
        args[pathkey] = os.path.abspath(args[pathkey])
    for listkey in ['--samples', '--patients', '--timepoints', '--years',
                    '--runs']:
        if args[listkey]:
            args[listkey] = [s.strip() for s in args[listkey].split(',')]


def selected(args, row):
    filters = [('--samples', 'SampleID'),
               ('--patients', 'PatientID'),
               ('--timepoints', 'TimePoint'),
               ('--years', 'SequencingYear'),
               ('--runs', 'RunID')]
    for optkey, colname in filters:
        if args[optkey] and row[colname] not in args[optkey]:
            return False
    return True


def main(args):
    if not os.path.exists(args['<outdir>']):
        os.makedirs(args['<outdir>'])
    os.chdir(args['<outdir>'])
    with open(args['<metadata_table>'], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        reader.fieldnames = [re.sub(r'^#\s*', '', name)
                             for name in reader.fieldnames]
        for row in reader:
            if not selected(args, row):
                continue

            sample_id = row['SampleID']
            demux_id = row['DemultiplexedID']
            logger.info(f"Processing sample {sample_id} "
                        f"(patient {row['PatientID']}, {row['TimePoint']})")
            logger.debug(f"Row info: {row}")

            if os.path.exists(sample_id) and not args['--force']:
                logger.warning(f"Output directory '{sample_id}' already "
                               "exists, skipping (use --force to override)")
                continue

            fastqsdir = (f"{args['<readsrootdir>']}/{row['RunID']}"
                         f"/{row['RunID']}_10X_RawData_Outs"
                         f"/{demux_id}/{row['FlowcellID']}")
            logger.debug(f"FASTQ directory: {fastqsdir}")

            cmdargs = ['count', '--id', sample_id, '--fastqs', fastqsdir,
                       '--transcriptome', args['<transcrdir>'],
                       '--localcores', args['--localcores'],
                       '--localvmem', args['--localvmem']]

            if args['--dry']:
                logger.info("Dry run")
                logger.info("Command: cellranger " + " ".join(cmdargs))
                cmdargs.append("--dry")

            sh.cellranger(*cmdargs, _fg=True)


args = docopt.docopt(__doc__, version='1.0')
postprocess_args(args)
setup_logger(args)
main(args)
