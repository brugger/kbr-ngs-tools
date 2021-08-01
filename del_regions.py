#!/usr/bin/env python3
""" 
 
 
 
 Kim Brugger (08 Jun 2021), contact: kim.brugger@uib.nok
"""

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

dels = {}

import pysam
samfile = pysam.AlignmentFile("tyt.sorted.bam", "rb" )
for pileupcolumn in samfile.pileup():
#    print ("\ncoverage at base %s = %s" %
#           (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        if pileupread.indel < -20 :
            read_name = pileupread.alignment.query_name
            if read_name not in dels:
                dels[ read_name ] = []
            dels[ read_name ].append( [pileupcolumn.pos+1, pileupcolumn.pos+1+ abs(pileupread.indel)] )
#            print( f"{pileupread.alignment.query_name} -- > {pileupread.indel}" )

pp.pprint( dels )
