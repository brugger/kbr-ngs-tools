#!/usr/bin/env python3

import os
import re
import sys
import argparse
from datetime import datetime, timedelta
import time
import json
import tempfile
import getpass
import pprint as pp

from pysam import VariantFile


def comp_vars( sample1:dict, sample2:dict) -> dict:

    total_vars = 0
    shared     = 0
    different  = 0
    missing    = 0

    for chrom in list(sample1.keys()):
        for pos in list(sample1[ chrom ].keys()):
            total_vars += 1
            if chrom not in sample2 or pos not in sample2[chrom]:
                missing += 1
            elif sample1[chrom][pos] != sample2[chrom][pos]:
                different += 1
                del sample2[chrom][pos]
            else:
                shared += 1
                del sample2[chrom][pos]
            

    for chrom in list(sample2.keys()):
        total_vars += 1
        for pos in list(sample2[ chrom ].keys()):
            if chrom not in sample1 or pos not in sample1[chrom]:
                missing += 1
            elif sample2[chrom][pos] != sample1[chrom][pos]:
                different += 1
            else:
                shared += 1


#    print(total_vars, shared, different, missing)

    if total_vars == 0:
        return 0.0

    return shared/total_vars*1.0
    


def main() -> None:

    samples = {}


    for infile in sys.argv[1:]:

        vcf_in = VariantFile(infile)  # auto-detect input format
        infile = re.sub(r'.*/', '', infile)
        infile = re.sub(r'\..*', '', infile)
#        print( infile )

        samples[ infile ] = {}

        for rec in vcf_in.fetch('chr22'):
            if len( rec.alts ) > 1:
                continue

            chrom = rec.chrom
            pos   = rec.pos
            ref   = rec.ref
            alt   = rec.alts[0]
            if chrom not in samples[ infile ]:
                samples[ infile ][ chrom ] = {}
            
            samples[ infile ][ chrom ][ pos ] = (ref, alt)

#            print( rec.chrom, rec.pos, rec.ref, alt )

    sample_names = list(samples.keys())
    for i in range(0, len(sample_names)-1):
        for j in range(i+1, len(sample_names)):
            shared_vars = comp_vars( samples[sample_names[ i ]], samples[sample_names[ j ]])
            print( f"{sample_names[ i ]}\t{sample_names[ j ]}\t{shared_vars}")




if __name__ == "__main__":
    main()

