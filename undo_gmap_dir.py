#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

parser.description = """
Undo GMAP's attempt at strand rectification in GMAP.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff )


#Will hold IDs of Dir=antisense genes & all their child features
antisense: set[str] = set()


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If this is a gene, check if it's antisense
	if entry.type_ == "gene":
		#The entry should have a "Dir" attribute
		assert "Dir" in entry.attrs
		
		#If this is an antisense gene, add it to the list
		if entry.attrs["Dir"] == "antisense":
			antisense.add( entry.attrs["ID"] )
	

	#If this is not a gene, check if it's a child of an antisense gene
	elif entry.attrs["Parent"] in antisense:
		antisense.add( entry.attrs["ID"] )
	
	
	#If needed, invert the strand & remove the Dir attribute
	if entry.attrs["ID"] in antisense:
		#'.' strand left untouched
		assert entry.strand != '.'
		if entry.strand == '+':
			entry.strand = '-'
		elif entry.strand == '-':
			entry.strand = '+'
	
	
	#Write entry
	print(entry)