#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr", metavar="ATTR",
					help = "Name of attribute to extract",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry possess the given attribute, report it
	if args.attr in entry.attrs:
		print(to_tsv( entry.attrs["ID"], entry.attrs[args.attr] ))
