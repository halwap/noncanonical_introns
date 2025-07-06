#!/usr/bin/env python3

from sys import argv
from argparse import ArgumentParser
from extras import *


parser = ArgumentParser()

#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
    help = "GFF file contatining annotation")

args = parser.parse_args()


require_files([args.gff])


#Loop over GFF records
for fields in read_tsv(args.gff, '#'):
	#Get gene records
	if fields[2] == "gene":
		
		#Convert attributes to dictionary
		attrs = attr_to_dict(fields[8])
		
		#Print gene ID & avg_intron_score if present
		print('\t'.join([ attrs["ID"], attrs.get("avg_intron_score", "") ]))
