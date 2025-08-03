#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *
from collections import defaultdict
from operator import attrgetter

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

#Flags
parser.add_argument("-0", "--zero",
					help = "Start counting from 0, not 1",
					action = "store_true", default = False)
parser.add_argument("-S", "--scaffold",
					help = "Order by position on scaffold, not within gene",
					action = "store_true", default = False)

parser.description = """
Generate new ordinal IDs for all 'exon' & 'intron' features.
"""

parser.epilog = """
Output can be passed to add_attrib.py to apply new IDs.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Will hold the child exons & introns of each mRNA feature
exons: dict[str,list[GFF]] = defaultdict(list)
introns: dict[str,list[GFF]] = defaultdict(list)


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If this is an exon or intron, submit it to the appropriate dictionary
	if entry.type_ == "exon":
		assert "Parent" in entry.attrs
		exons[entry.attrs["Parent"]].append( entry )
	
	if entry.type_ == "intron":
		assert "Parent" in entry.attrs
		introns[entry.attrs["Parent"]].append( entry )


#For each mRNA, sort its exons appropriately & generate new IDs
for parent, children in exons.items():
	#Sort according to position on scaffold if explicitly requested, or if this is a positive
	#strand mRNA
	if args.scaffold or children[0].strand == '+':
		children.sort(key=attrgetter("start"))
	#Sort opposite to position on scaffold is this is a negative strand mRNA, and within-gene
	#ordering was requested
	else:
		children.sort(key=attrgetter("start"), reverse=True)
	
	#Generate & print new IDs
	for idx, exon in enumerate(children, 0 if args.zero else 1):
		print(to_tsv( exon.attrs["ID"], f"{parent}.exon{idx}" ))


#Likewise for introns
for parent, children in introns.items():
	if args.scaffold or children[0].strand == '+':
		children.sort(key=attrgetter("start"))
	else:
		children.sort(key=attrgetter("start"), reverse=True)
	
	for idx, intron in enumerate(children, 0 if args.zero else 1):
		print(to_tsv( intron.attrs["ID"], f"{parent}.intron{idx}" ))
