#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("type", metavar="TYPE",
					help = "New feature type",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "List of features to modify",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


#Get new value of attribute for listed features
list_: set[str] = set( parse_list( args.list ) )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry is a listed one, set its feature type
	if entry.attrs["ID"] in list_:
		entry.type_ = args.type
	print( entry )
