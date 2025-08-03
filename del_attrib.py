#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr", metavar="ATTR",
					help = "Name of attribute to delete",
					type = str)
parser.add_argument("list", metavar="LIST", nargs='?',
					help = "List of features to modify",
					type = str)

parser.description = """
Process GFF to remove attribute ATTR based on LIST.
"""

parser.epilog = """
LIST is a listfile of feature IDs.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


if args.list:
	ids: set[str] = set( parse_list( args.list ) )


for entry in parse_gff(args.gff):
	if ( args.list is None or entry.attrs["ID"] in ids ) and args.attr in entry.attrs:
		del entry.attrs[args.attr]
	print( entry )
