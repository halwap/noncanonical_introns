#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *
from extras import merge_sets

###############################################################################

parser = ArgumentParser()

#Positional argument
parser.add_argument("pairs", metavar="PAIRS",
					help = "List of pairs of features",
					type = str)
parser.add_argument("lengths", metavar="LENGTHS",
					help = "List of feature lengths",
					type = str)

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
