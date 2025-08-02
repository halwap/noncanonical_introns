#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr", metavar="ATTR", nargs='?',
					help = "Attribute to optionally report",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If a specific attribute was requested, and the given attribute posesses it
	if args.attr and args.attr in entry.attrs:
		#Report ID & requested attribute
		print(to_tsv( entry.attrs["ID"], entry.attrs[args.attr] ))
	#If no attribute was requested, report just the ID
	elif not args.attr:
		print(entry.attrs["ID"])