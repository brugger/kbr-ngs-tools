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

from pysam import VariantFile, TabixFile


def main() -> None:


    parser = argparse.ArgumentParser(description=f'classify probes into good and bad (true/false)')

    parser.add_argument('-d', '--depth', help="min variant depth ", default=10)
    parser.add_argument('-q', '--quality', help="min quality score ", default=50)
    parser.add_argument('-c', '--chrom', help="run single chromosome, mainly for development ", default=None)
    parser.add_argument('-g', '--gq', help="sample genotype quality ", default=0)
    parser.add_argument('-s', '--snps-only', help="Only consider SNPS ", default=True, action='store_false')
    parser.add_argument('-m', '--multi-vars', help="add multi variant sites to analysis  ", default=False, action='store_true')
    parser.add_argument('-r', '--reps-file', help='tab file with rep regions (mq=0)')
    parser.add_argument('infile', nargs=1, help="file to analyse")

    args = parser.parse_args()

    args.infile = args.infile[0]

    total_vars = 0
    low_depth  = 0
    low_quality = 0
    concordant = 0
    discordant = 0
    rep_region = 0
    no_call    = 0
    low_gq     = 0
    odd_homozygosity = 0


    tbx = None
    if args.reps_file:
        tbx = TabixFile(args.reps_file)

    vcf_in = VariantFile(args.infile)  # auto-detect input format
    for rec in vcf_in.fetch(args.chrom):

        if  len( rec.alts ) > 1 and args.multi_vars == False:
            continue

        if  len( rec.alts[0]) > 1:
            continue

        if rec.qual < int( args.quality ):
            low_quality += 1
            continue

        total_vars += 1

        min_depth = -1
        for dp in [rec.samples[s]['DP'] for s in rec.samples ]:
            if dp is None:
                dp = 0
            if min_depth == -1 or min_depth > dp:
                min_depth = dp

        if min_depth < int( args.depth ):
#            print(f" Low: {min_depth}")
            low_depth += 1
            continue

        gq_low = False
        for gq in [rec.samples[s]['GQ'] for s in rec.samples ]:
            if gq is not None and gq < int(args.gq):
                low_gq += 1
                gq_low = True
                continue

        if gq_low:
            continue

        if tbx is not None:
            rows = tbx.fetch(rec.chrom, rec.pos, rec.pos+1)
            # Magic that checks if the iterator is empty or not!
            if any( True for _ in rows ):
                rep_region += 1
                continue


        if rec.samples[0]['GT'] == rec.samples[1]['GT']:
            concordant += 1
        elif None in rec.samples[0]['GT'] + rec.samples[1]['GT']:
            no_call += 1
        elif rec.samples[0]['GT'] == (1,1) and  rec.samples[1]['GT'] == (0,1):
            print( f"{rec.chrom}:{rec.pos}\t{rec.samples[0]['GT']} != {rec.samples[1]['GT']}" )
            odd_homozygosity +=1
            discordant += 1

        else:
            discordant += 1


    print(f"total: {total_vars}, concordant: {concordant}, discordant: {discordant}")
    print(f"Filtered out: low-depth: {low_depth}, no_call: {no_call}, rep_regions: {rep_region}, low_gq: {low_gq}, low quality: {low_quality}, odd_homozygosity: {odd_homozygosity}")
    print(f"corcondance: {(concordant*100)/(discordant+concordant)} %")


if __name__ == "__main__":
    main()

