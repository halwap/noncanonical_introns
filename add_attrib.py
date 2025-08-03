#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr", metavar="ATTR",
					help = "Name of attribute to add",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "List of features to modify, w/ new attribute values",
					type = str)

parser.description = """
Process GFF to add/modify attribute ATTR based on LIST.
"""

parser.epilog = """
LIST is a 2-column TSV. Field 1 is a feature's ID, field 2 is the new value
of ATTR to write for said feature.

Any path can be '-' to read from stdin. Writes to stdout.
"""

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
