#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *
from extras import merge_sets

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional argument
parser.add_argument("pairs", metavar="PAIRS",
					help = "List of pairs of features",
					type = str)
parser.add_argument("lengths", metavar="LENGTHS",
					help = "List of feature lengths",
					type = str)

parser.description = """
Given a list of relations between features in PAIRS, and the length of each
feature in LENGTHS, find and report each feature's most representative related
feature, if it has one.
"""

parser.epilog = """
PAIRS is a 2-column TSV, both columns containing feature IDs. If two features
are listed on the same line, this signifies they are related.
This can be the output of find_altseq_gff.py.

LENGTHS is a 2-column TSV. Field 1 is a feature's ID, field 2 is its length.
This can be the output of get_ft_len.py.

Output can be passed to apply_adoptions.py.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.pairs, args.lengths )


#Get length of each features
lengths: dict[str,int] = dict( parse_tsv( args.lengths ) )

#Get a list of features pairings
pairs: list[set[str]] = [ set( ids ) for ids in parse_tsv( args.pairs ) ]

#Merge pairs into complete groups
groups: list[set[str]] = merge_sets(pairs)


#Iterate over groups
for group in groups:
	
	#Get the longest feature in each group
	longest: str = max(group, key = lambda x: lengths[x])
	

	#For each feature in the group, report that it's children should be adopted by `longest',
	#except for `longest' itself
	for id_ in group:
		if id_ != longest:
			print(to_tsv( id_, longest ))
