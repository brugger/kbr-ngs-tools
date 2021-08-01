#!/usr/bin/env python3
""" 
 
 
 
 Kim Brugger (08 Jun 2021), contact: kim.brugger@uib.nok
"""

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

import pysam

BAM_CMATCH	= 0
BAM_CINS	= 1
BAM_CDEL	= 2
BAM_CREF_SKIP	= 3
BAM_CSOFT_CLIP	= 4
BAM_CHARD_CLIP	= 5
BAM_CPAD	= 6
BAM_CEQUAL	= 7
BAM_CDIFF	= 8
BAM_CBACK	= 9



samfile = pysam.AlignmentFile("tyt.sorted.bam", "rb")
print( samfile.header, end='')

event = ['', 'I', 'D']

for read in samfile.fetch():
     if not read.is_unmapped:
         cigartubles = read.cigartuples
         for ct in cigartubles:
              if ct[0] in [BAM_CDEL]  and ct[1] >= 50:
                   print( read.to_string() )
                   break
#                   print( f"{event[ct[0]]}\t{ct[1]}\t{read.query_name}");
#                   print( f">{read.query_name}")
#                   print( read.query_sequence )
