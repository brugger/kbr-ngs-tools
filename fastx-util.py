#!/usr/bin/env python3
""" 
 
 
 
 Kim Brugger (08 Jun 2021), contact: kim.brugger@uib.no
"""

import sys
import re
import argparse
import statistics
import pprint
pp = pprint.PrettyPrinter(indent=4)

import pysam

import kbr.args_utils as args_utils
import kbr.version_utils as version_utils
import kbr.string_utils as string_utils




def fastx_grep(args):

    filename = args_utils.get_or_fail(args, "filename to use is required")
    pattern = args_utils.get_or_fail(args, "filename to use is required")
    pattern = re.compile( pattern )
    with pysam.FastxFile(filename) as fh:
        for entry in fh:
            if pattern.search( entry.name ) or pattern.search(entry.sequence):
                print( entry )


def fastx_stats(args, lengths_only:bool=False, counts_only:bool=False):

    filename = args_utils.get_or_fail(args, "filename to use is required")

    count = 0
    lengths = []
    
    with pysam.FastxFile(filename) as fh:
        for entry in fh:
            count += 1
            if counts_only:
                continue
        
            if lengths_only:
                print(f"{entry.name}\t{len(entry.sequence)}")
                continue

            lengths.append( len(entry.sequence) )
    
    if counts_only:
        print(f"{count} entries")
    if lengths_only:
        pass
    else:
        lengths = sorted( lengths )
        print(f"Mean length: {statistics.mean(lengths):.2f}")
        print(f"Median length: {statistics.median(lengths):.2f}")
        print(f"Mode length: {statistics.mode(lengths):.2f}")
        print(f"Std dev: {statistics.stdev(lengths):.2f}")
        print(f"Shortest: {lengths[0]}")
        print(f"Longest: {lengths[-1]}")
        print("quantiles")
        quantiles = statistics.quantiles(lengths, n=20)
        for i in range(0, len( quantiles )):
            print( f"{(i+1)*5}%\t{quantiles[i]}")
                       
        
def main():


    commands = {'g': 'grep', 'l': 'lengths', 'c': 'count', 's': 'stats', 'h':'help'}
    parser = argparse.ArgumentParser(description=f'fastx-utils: fasta/q utils')

    parser.add_argument('command', nargs='*', help="{}".format(",".join(commands)))

    args = parser.parse_args()

    args_utils.min_count(1, len(args.command),
                         msg="fastx-utils takes one of the following commands: {}".format(string_utils.comma_sep(list(commands.values()))))

    command = args.command.pop(0)


    if command in commands:
        command = commands[ command ]

    if command == 'grep':
        fastx_grep(args.command)
    elif command == 'lengths':
        fastx_stats(args.command, lengths_only=True)
    elif command == 'count':
        fastx_stats(args.command, counts_only=True)
    elif command == 'stats':
        fastx_stats(args.command)
    else:
        print("The tool support the following commands: {}\n".format(string_utils.comma_sep(list(commands.values()))))
        parser.print_usage()
        parser.add_argument('command', nargs='+', help="{}".format(",".join(commands)))
        sys.exit(1)

if __name__ == "__main__":
    main()

        
