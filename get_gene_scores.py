#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

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
	#Get gene records
	if entry.type_ == "gene":
		#Print gene ID & avg_intron_score if present
		print(to_tsv( entry.attrs["ID"], entry.attrs.get("avg_intron_score", "NA") ))
