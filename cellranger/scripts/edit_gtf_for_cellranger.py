#!/usr/bin/env python3

import sys
import re

with open(sys.argv[1]) as f:
    for orig_line in f:
        if orig_line.startswith('#'):
            continue
        line = orig_line.strip().split("\t")
        feature = line[2]
        if feature == 'gene':
            print(orig_line.strip())
            chrom = line[0]
            source = line[1]
            start = line[3]
            end = line[4]
            score = line[5]
            strand = line[6]
            frame = line[7]
            attributes = [attr.strip().split(" ", 1) for attr in line[8].split(";")]
            gene_id = [attr[1].strip('"') \
                    for attr in attributes if attr[0] == 'gene_id'][0]
            transcript_id = 't_' + gene_id
            has_trascript_id = [attr[0] for attr in attributes \
                    if attr[0] == 'transcript_id']
            transcript_id_attr = ['transcript_id', '"'+transcript_id+'"']
            if has_trascript_id:
                attributes = [attr if attr[0] != 'transcript_id' \
                        else transcript_id_attr for attr in attributes]
            else:
                attributes.append(transcript_id_attr)
            print("\t".join([chrom, source, 'transcript', start, end, score, \
                    strand, frame, '; '.join([' '.join(attr) \
                    for attr in attributes])]))
            print("\t".join([chrom, source, 'exon', start, end, score, \
                    strand, frame, '; '.join([' '.join(attr) \
                    for attr in attributes])]))
