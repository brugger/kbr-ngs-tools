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

import pysam

def pseudo_region(quals:list) -> bool:
    for qual in quals:
        if qual == 0:
            return True

    return False

def avg_qual(quals:list) -> int:
    summed = 0
    if len(quals) == 0:
        return -1
    for qual in quals:
        summed += qual


    return summed/len(quals)




def main() -> None:

    
    samfile = pysam.AlignmentFile("NA12878_FN.bam", "rb")

    SNPs = []

    file_handle = open("PP_FN_NO_ALTS.entries", 'r')
    for l in  file_handle.readlines():

        chrom, pos, ref, alt = l.rstrip().split("\t")
        pos = int(pos)
        pseudo = False
        var_type = f"SNP\t{ref}/{alt}"
        if len(ref) > 1 or len(alt) > 1:
            var_type = f"INDEL\t{ref}/{alt}"

        for pileupcolumn in samfile.pileup(chrom, pos - 1, pos):
            if pseudo_region ( pileupcolumn.get_mapping_qualities()):
                print(f"{var_type}\t{var_type}\t{chrom}:{pos}\tcoverage {pileupcolumn.n}\tPseudo region")
                pseudo = True
                break

        if not pseudo:
            for pileupcolumn in samfile.pileup(chrom, pos - 1, pos):
                if pileupcolumn.pos != pos -1:
                    continue
                bases = {'A':0,'C':0,'G':0,'T':0,'N':0}
            # query position is None if is_del or is_refskip is set.

                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        bases[pileupread.alignment.query_sequence[pileupread.query_position]] += 1


                avg_mq = avg_qual(pileupcolumn.get_mapping_qualities())
                if pileupcolumn.n < 10:
                    print(f"{var_type}\t{chrom}:{pos}\tcoverage {pileupcolumn.n}\tavg mq: {avg_mq:.2f}\tlow coverage")
                    pseudo = True
                else:
                    print(f"{var_type}\t{chrom}:{pos}\tcoverage {pileupcolumn.n}\tavg mq: {avg_mq:.2f}\tA:{bases['A']}/C:{bases['C']}/G:{bases['G']}/T:{bases['T']}/N:{bases['N']}")
                    pseudo = True
        
        if not pseudo:
            print(f"{var_type}\t{chrom}:{pos} no coverage")

if __name__ == "__main__":
    main()

