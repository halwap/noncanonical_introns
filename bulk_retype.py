#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "List of features to modify, w/ new type name",
					type = str)

parser.description = """
Process GFF to modify the type of features as per LIST.
"""

parser.epilog = """
LIST is a 2-column TSV. Field 1 is a feature's ID, field 2 is the new name
of the feature's type.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


#Get new type name for listed features
types: dict[str,str] = dict( parse_tsv( args.list ) )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry is a listed one, modify its type
	if entry.attrs["ID"] in types:
		entry.type_ = types[entry.attrs["ID"]]
	print( entry )
