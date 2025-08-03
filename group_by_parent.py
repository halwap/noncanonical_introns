#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *
from collections import defaultdict

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

parser.description = """
Group features in GFF, such that two features end up in the same group if they
have the same parent.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Will hold the children of each parent
children: dict[str:list[str]] = defaultdict(list)


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If the entry is a child, add it to the list of its Parents
	if "Parent" in entry.attrs:
		children[entry.attrs["Parent"]].append( entry.attrs["ID"] )


#Print each parent's children on separate lines
for child_list in children.values():
	print(to_tsv( *child_list ))