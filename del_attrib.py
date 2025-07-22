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
					help = "File containing IDs to modify",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


#Get new value of attribute for listed features
attrs: set[str] = set( parse_list( args.list ) )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry is a listed one and has the given attribute, remove it
	if entry.attrs["ID"] in attrs and args.attr in entry.attrs:
		del entry.attrs[args.attr]
	print( entry )
