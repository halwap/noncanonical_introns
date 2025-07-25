#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr1", metavar="ATTR1",
					help = "Name of attribute to modify",
					type = str)
parser.add_argument("attr2", metavar="ATTR2",
					help = "New name for the attribute",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry has the appropriate argument, rename it
	if args.attr1 in entry.attrs:
		entry.attrs[args.attr2] = entry.attrs[args.attr1]
		del entry.attrs[args.attr1]
	print( entry )
