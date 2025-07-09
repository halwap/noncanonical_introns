#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *
from libintrons import CONV_SS

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")


args = parser.parse_args()



###############################################################################
require_files( args.gff )


#Loop over GFF records
for entry in parse_gff(args.gff):
	#If this record is an intron
	if entry.type_ == "intron":
		
		if "splice_site" in entry.attrs:
			#Add _C or _N depending on splice site
			entry.type_ += "_C" if entry.attrs["splice_site"][:4] in CONV_SS else "_N"
	
	print(entry)
