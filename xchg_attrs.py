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
					help = "Name of one attribute",
					type = str)
parser.add_argument("attr2", metavar="ATTR2",
					help = "Name of another attribute",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#Get the appropriate args
	attr1 = entry.attrs.get(args.attr1, None)
	attr2 = entry.attrs.get(args.attr2, None)
	
	#Swap
	if attr1 is not None and attr2 is not None:
		entry.attrs[attr1] = attr2
		entry.attrs[attr2] = attr1
	
	elif attr1 is not None:
		entry.attrs[attr2] = attr1
		del entry.attrs[attr1]

	elif attr2 is not None:
		entry.attrs[attr1] = attr2
		del entry.attrs[attr2]

	print( entry )
