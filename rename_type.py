#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("prev", metavar="TYPE1",
					help = "Old feature type",
					type = str)
parser.add_argument("new", metavar="TYPE2",
					help = "New feature type",
					type = str)

parser.description = """
Rename type TYPE1 to TYPE2 in GFF.
"""

parser.epilog = """
Any path can be '-' to read from stdin or write to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


for entry in parse_gff(args.gff):
	if entry.type_ == args.prev:
		entry.type_ = args.new
	print( entry )
