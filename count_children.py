#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *
from collections import defaultdict

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

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
