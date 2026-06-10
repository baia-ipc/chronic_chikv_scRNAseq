#!/usr/bin/env python3
"""
Collect metric summaries for all samples in the output directories of cellranger count

Usage:
  ./collect_metric_summaries.py [options] <outputdir>

Arguments:
  <outputdir>  Directory containing output directories of cellranger count

Options:
    -m, --metadata FILE  Metadata file
    -c, --colnames TSV   Remapping column names before output
                         (TSV file oldname<tab>newname)
    -v, --verbose        Print more output
    -d, --debug          Print debug output
    -q, --quiet          Print less output
    -h, --help           Show this screen
    -V, --version        Show version
"""
import csv
import docopt
import sys
from pathlib import Path
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

def parse_metrics_summary(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)
        values = next(reader)
        for i in range(len(values)):
            if values[i].endswith('%'):
                values[i] = float(values[i][:-1])/100
            elif ',' in values[i]:
                values[i] = int(values[i].replace(',',''))
            else:
                values[i] = int(values[i])
        return dict(zip(header, values))

def main(args):
    # parse all metrics_summary.csv files
    # contained under <outputdir>/<sample>/outs
    metrics = {}
    for sample in Path(args['<outputdir>']).glob('*'):
        if sample.is_dir():
            metrics[sample.name] = \
                parse_metrics_summary(sample / 'outs' / 'metrics_summary.csv')
    # parts metadata file
    if args['--metadata']:
        with open(args['--metadata'], 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            # get the name of the leftmost column
            sample_col = reader.fieldnames[0]
            for row in reader:
                # remove sample_col from row
                sample_id = row.pop(sample_col)
                # merge with metric[sample] if sample is in metrics
                if sample_id in metrics.keys():
                    row.update(metrics[sample_id])
                    metrics[sample_id] = row
    # parse colnames file with lines <oldname>\t<newname>
    if args['--colnames']:
        with open(args['--colnames'], 'r') as f:
            colnames = {}
            for line in f:
                oldname, newname = line.strip().split('\t')
                colnames[oldname] = newname
    # header
    print('sample', end='')
    for key in metrics[sample.name].keys():
        if args['--colnames'] and key in colnames.keys():
            key = colnames[key]
        print('\t'+key, end='')
    print()
    # values
    for sample in metrics.keys():
        print(sample, end='')
        for key in metrics[sample].keys():
            print('\t'+str(metrics[sample][key]), end='')
        print()

args = docopt.docopt(__doc__, version='0.1')
setup_logger(args)
main(args)
