#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("prev", metavar="TYPE1",
					help = "Old feature type",
					type = str)
parser.add_argument("new", metavar="TYPE1",
					help = "New feature type",
					type = str)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


for entry in parse_gff(args.gff):
	if entry.type_ == args.prev:
		entry.type_ = args.new
	print( entry )
