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



samfile = pysam.AlignmentFile("test_reads_hangs.bam", "rb")

for read in samfile.fetch():
     if not read.is_unmapped:
         cigar = read.cigartuples
         if read.query_sequence is None:
             continue
         if cigar[0][0] == BAM_CSOFT_CLIP and cigar[0][1] > 10:
             print( f">{read.query_name}-us")
             print( read.query_sequence[0: cigar[0][1]] )
             if read.query_qualities:
                 print( read.query_qualities[0: cigar[0][1]] )
        
         if cigar[-1][0] == BAM_CSOFT_CLIP and cigar[-1][1] > 10:
             print( f">{read.query_name}-ds")
             print( read.query_sequence[-cigar[-1][1]: -1] )
             if read.query_qualities:
                 print( read.query_qualities[-cigar[-1][1]: -1] )
             
#             sys.exit()
