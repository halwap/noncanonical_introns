#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("list", metavar="LIST", nargs='?',
					help = "List of types to report",
					type = str)

parser.description = """
List all features in GFF of one of given types.
"""

parser.epilog = """
LIST is a listfile of feature types.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


if args.list:
	types: set[str] = set( parse_list( args.list ) )


for entry in parse_gff(args.gff):
	if args.list is None or entry.type_ in types:
		print( entry.attrs["ID"] )
