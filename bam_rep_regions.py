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

import pysam



def main() -> None:

    infile = sys.argv[1]
    samfile = pysam.AlignmentFile(infile, "rb" )
    region_start = None
    last_pos     = None
    last_chrom   = None
    for pileupcolumn in samfile.pileup("chr1"):

        if region_start is not None and (last_chrom != pileupcolumn.reference_name or pileupcolumn.reference_pos != last_pos + 1):
            if last_pos is not None:
                print("\t".join( [last_chrom, str(region_start), str(last_pos)] ))
            region_start = None

        if 0 in pileupcolumn.get_mapping_qualities():
            if region_start is None:
                region_start = pileupcolumn.reference_pos
                last_chrom = pileupcolumn.reference_name

            last_pos = pileupcolumn.reference_pos

    print("\t".join( [last_chrom, str(region_start), str(last_pos)] ))

    samfile.close()




if __name__ == "__main__":
    main()
