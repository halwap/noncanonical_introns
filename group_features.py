#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *
from extras import merge_sets
from itertools import combinations

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

#Flags
parser.add_argument("-t", "--type", metavar="TYPE",
					help = "Only consider features of type TYPE",
					type = str, default = "")
parser.add_argument("-S", "--stranded",
					help = "Include strand information when grouping",
					action = "store_true", default = False)
parser.add_argument("-r", "--rough",
					help = "Group features based only on a shared seqid and optionally strand",
					action = "store_true", default = False)
parser.add_argument("-M", "--multi-elem",
					help = "Only print groups containing multiple elements",
					action = "store_true", default = False)

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Get GFF entries
if args.type:
	entries = list( entry for entry in parse_gff(args.gff) if entry.type_ == args.type )
else:
	entries = list( entry for entry in parse_gff(args.gff) )

#If no matching entries were found, exit early
if not entries:
	eprint("No matching entries found in GFF file")
	exit()


#Will hold every pair of 2 entries which should end up in the same group
pairs: list[tuple[str,str]] = []


#Iterate over all unique pairs of entries in `entries'
for entry1,entry2 in combinations(entries, 2):
	#Determine if the two entries should belong in the same group
	same_group = False
	if args.rough:
		if entry1.seqid == entry2.seqid and (not args.stranded or entry1.strand == entry2.strand ):
			same_group = True
	else:
		same_group = entry1.overlaps(entry2, unstranded = not args.stranded)
	
	
	#If the two entries should belong in the same group, submit them to `pairs'
	if same_group:
		pairs.append( (entry1.attrs["ID"], entry2.attrs["ID"]) )


#Merge pairs into complete groups
groups: list[set[str]] = merge_sets(pairs)


#Report groups
for group in groups:
	print(to_tsv( *group ))


#If single-element groups were requested, print them also
if not args.multi_elem:
	#Get the set of all IDs
	ids: set[str] = set( entry.attrs["ID"] for entry in entries )
	
	#Remove IDs already printed
	for group in groups:
		ids -= group
	
	
	#Print all remaining IDs as remaining single-element groups
	for id_ in ids:
		print(id_)
