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
Count the number of children each feature in GFF has.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Will hold the number of children of each feature
child_cnt: dict[str,int] = defaultdict(int)


#Loop over GFF records
for entry in parse_gff(args.gff):
	#Instantiate this ID's record if necessary
	child_cnt[entry.attrs["ID"]]
	
	#Count the current feature, if it is a child
	if "Parent" in entry.attrs:
		child_cnt[entry.attrs["Parent"]] += 1


#Print children counts
for id_, cnt in child_cnt.items():
	print(to_tsv( id_, cnt ))
