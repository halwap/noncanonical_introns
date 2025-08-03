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
Get length of every feature in GFF.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


lengths: dict[str,int] = defaultdict(int)


#Loop over GFF records
for entry in parse_gff(args.gff):
	lengths[entry.attrs["ID"]] += len(entry)

#Report each feature's length.
for id_, length in lengths.items():
	print(to_tsv( id_, length ))
