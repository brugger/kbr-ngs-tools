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
import copy
import pprint as pp

from pysam import VariantFile


def read_args(filename:str) -> list:

    if filename == '-':
        file_handle = sys.stdin
    else:
        file_handle = open(filename, 'r')

    content = file_handle.read()
    file_handle.close()

    args = re.sub("\n", ' ', content).split()
    return args




def comp_vars( sample1:dict, sample2:dict) -> int:

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


#    print(f"total: {total_vars}, shared: {shared}, differenct: {different}, missing: {missing}")

    if total_vars == 0:
        return 0

    return shared, different, missing



def readin_samples( infiles:list, chrom:str='' )-> dict:

    samples = {}

    for infile in infiles:

        vcf_in = VariantFile(infile)  # auto-detect input format
        infile = re.sub(r'.*/', '', infile)
        infile = re.sub(r'\..*', '', infile)

        samples[ infile ] = {}

        for rec in vcf_in.fetch(chrom):
            if len( rec.alts ) > 1:
                continue

            chrom = rec.chrom
            pos   = rec.pos
            ref   = rec.ref
            alt   = rec.alts[0]
            if chrom not in samples[ infile ]:
                samples[ infile ][ chrom ] = {}

            samples[ infile ][ chrom ][ pos ] = (ref, alt)


    return samples



def main() -> None:


    parser = argparse.ArgumentParser(description=f'intersect loads on vcfs files')

    parser.add_argument('-r', '--reference-samples', help="reference samples ")
    parser.add_argument('-s', '--samples-file', help="samples ")
    parser.add_argument('-c', '--chrom', help="samples ", default='chr22')

    parser.add_argument('samples', nargs='*', help="samples to compare against reference samples")

    args = parser.parse_args()

    args.reference_samples =  read_args( args.reference_samples)
    if args.samples_file is not None:
        args.samples +=  read_args( args.samples_file)

    print( 'reading in refs')
    refs    = readin_samples( args.reference_samples, args.chrom)

    print( 'analysing samples')
    for sample in args.samples:
        vars = readin_samples( [sample], args.chrom)
        vars = vars[ list(vars.keys())[0]]

        score = 0
        name  = None
        for ref in refs:
            shared_vars, different, missing = comp_vars( copy.deepcopy(vars), copy.deepcopy(refs[ref]))
            if shared_vars > score:
                score = shared_vars
                name = ref

        print( f"{sample}\t{name}\t{score}")


if __name__ == "__main__":
    main()

