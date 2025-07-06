#!/usr/bin/env python3

from sys import argv
from argparse import ArgumentParser
from extras import *

###############################################################################
parser = ArgumentParser()

#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

args = parser.parse_args()



###############################################################################
require_files([args.gff])


#Loop over GFF records
for fields in read_tsv(args.gff, '#'):
	#If this record is an intron
	if fields[2] == "intron":
		
		#Convert attributes to dictionary
		attrs = attr_to_dict(fields[8])
		
		if "splice_site" in attrs:
			splice_site = attrs["splice_site"][:4]
			
			#Add _C or _N depending on splice site
			fields[2] += "_C" if splice_site in { "GTAG", "GCAG", "CTAC", "CTGC" } else "_N"
	
	print('\t'.join(fields))
