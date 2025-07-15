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
					help = "Name of attribute to add",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "File containing IDs to modify with new attribute values",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


#Get new value of attribute for listed features
attrs: dict[str,str] = dict( parse_tsv( args.list ) )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry is a listed one, modify its attributes
	if entry.attrs["ID"] in attrs:
		entry.attrs[args.attr] = attrs[entry.attrs["ID"]]
	print( entry )
