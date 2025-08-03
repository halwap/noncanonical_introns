#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

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

parser.description = """
Modify GFF to change the types of all features in LIST to TYPE.
"""

parser.epilog = """
LIST is a listfile of feature IDs.

Any path can be '-' to read from stdin. Writes to stdout.
"""

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
