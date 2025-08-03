#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("attr", metavar="ATTR", nargs='?',
					help = "Attribute to optionally report",
					type = str)

parser.description = """
Report value of attribute ATTR for every feature in GFF which has it.
"""

parser.epilog = """
If ATTR is not passed, all feature IDs are listed instead.

Any path can be '-' to read from stdin. Writes to stdout.
"""

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