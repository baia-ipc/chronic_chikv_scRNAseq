#!/usr/bin/env python3
"""
Analyze a SAM file and output UMI counts and strand information for each cell
and feature. With the --consolidate switch, it provides a compact report.

Usage:
    strand_analysis.py [options] <filename> [<feature>...]

Arguments:
    filename  The path to the SAM or BAM file.
    feature   The features to analyze. Multiple features can be specified.
              (default: all)

Options:
    --help                  Show this message.
    --version               Show version information.
    --consolidate           Output a compact summary report.
    -r, --region <region>   Only analyze reads in the specified region.
                            For performance reasons is recommended instead:
                            samtools view (-b <bam>|<sam>) <region>
"""

import pysam
import sys
from docopt import docopt
from collections import defaultdict, Counter

def get_strand_status(strand_info, mixed_simple = False):
    pos = strand_info['positive']
    neg = strand_info['negative']
    if pos > 0 and neg == 0:
        return "positive"
    elif neg > 0 and pos == 0:
        return "negative"
    elif mixed_simple:
        return "mixed"
    else:
        return "mixed (pos: {}, neg: {})".format(pos, neg)

def detailed_report(bfus_counts):
    for barcode, fus_counts in bfus_counts.items():
        print(f"Barcode: {barcode}")
        for feature, us_counts in fus_counts.items():
            print(f"  Feature: {feature}")
            for umi, s_counts in us_counts.items():
                count = s_counts['positive'] + s_counts['negative']
                strand_status = get_strand_status(s_counts)
                print(f"    UMI: {umi}, Count: {count}, Strand: {strand_status}")

def consolidated_report(bfus_counts):
    n_cells = len(bfus_counts)
    count_barcodes_by_feature_set = defaultdict(int)
    count_umis_by_feature_set = defaultdict(list)
    count_umis_by_strand_status = defaultdict(lambda: \
                                              defaultdict(lambda: \
                                                          defaultdict(int)))
    for barcode, fus_counts in bfus_counts.items():
        feature_set = set()
        for feature, us_counts in fus_counts.items():
            feature_set.add(feature)
        feature_set_name = "{} only".format(feature) if len(feature_set) == 1 \
                            else " and ".join(feature_set)
        count_barcodes_by_feature_set[feature_set_name] += 1 
        for feature, us_counts in fus_counts.items():
            umi_count = 0
            umi_count_by_strand_status = defaultdict(int)
            for umi, s_counts in us_counts.items():
                umi_count += 1
                strand_status = get_strand_status(s_counts, mixed_simple=True)
                umi_count_by_strand_status[strand_status] += 1
            count_umis_by_feature_set[feature_set_name].append(umi_count)
            count_umis_by_strand_status[feature_set_name][feature]['positive'] +=\
                  umi_count_by_strand_status['positive']
            count_umis_by_strand_status[feature_set_name][feature]['negative'] +=\
                  umi_count_by_strand_status['negative']
            count_umis_by_strand_status[feature_set_name][feature]['mixed'] +=\
                  umi_count_by_strand_status['mixed']
    print(f"Cells: {n_cells}")
    feature_set_names = sorted(count_barcodes_by_feature_set.keys(), 
                                key=lambda x: (x.endswith(' only'), x))
    for feature_set_name in feature_set_names:
        is_only = feature_set_name.endswith(' only')
        counter = Counter(count_umis_by_feature_set[feature_set_name])
        counter_out = []
        for c, n in counter.items():
            counter_out.append(f"{c} cells {n} UMIs")
        counter_str = ", ".join(counter_out)
        print(f"  {feature_set_name}: "+\
              f"{count_barcodes_by_feature_set[feature_set_name]} cells; "+\
              f"{sum(count_umis_by_feature_set[feature_set_name])} UMIs" +\
              f" [{counter_str}]", end='')
        for feature, strand_counts in \
                count_umis_by_strand_status[feature_set_name].items():
            print(" -- ", end='')
            if not is_only:
                print(f"{feature}: ", end='')
            to_print = []
            for strand_status, short_status in \
                    [('positive', 'pos'), ('negative', 'neg'), ('mixed', 'mix')]:
                if strand_counts[strand_status] > 0:
                    to_print.append(f"{strand_counts[strand_status]} {short_status}")
            print(", ".join(to_print), end='')
        print()

def main(file_path, features=None, region=None, consolidate=False):
    bfus_counts = \
       defaultdict(lambda: \
          defaultdict(lambda: \
            defaultdict(lambda: Counter({'positive': 0, 'negative': 0}))))

    with pysam.AlignmentFile(file_path, "r", check_sq=False) as samfile:
        for read in samfile:
            if region and not read.reference_name == region:
                continue
            feature = read.get_tag("GX") \
                        if read.has_tag("GX") else "Not-Assigned-to-Feature"
            #if feature is None:
            #    print("Warning: aln without GX tag", file=sys.stderr)
            #    print("  ", end='', file=sys.stderr)
            #    print(read.to_string(), file=sys.stderr)
            #    continue
            if features and feature not in features:
                continue
            if read.has_tag("UB") and read.has_tag("CB"):
                umi = read.get_tag("UB")
                barcode = read.get_tag("CB")
                strand = 'positive' if not read.is_reverse else 'negative'
                bfus_counts[barcode][feature][umi][strand] += 1

    if consolidate:
        consolidated_report(bfus_counts)
    else:
        detailed_report(bfus_counts)

if __name__ == "__main__":
    arguments = docopt(__doc__, version="1.0")
    main(arguments["<filename>"], arguments["<feature>"],
         arguments["--region"], arguments["--consolidate"])

