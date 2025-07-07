#!/usr/bin/env python3

from argparse import ArgumentParser
from extras import *

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

parser.add_argument("list", metavar="LIST",
	help = "File containing features to discard")


#Options
parser.add_argument("-k", "--keep",
	help = "Keep listed features, and discard unlisted ones",
	default = False, action = "store_true")


args = parser.parse_args()



###############################################################################
require_files( args.gff, args.list )


#Deserialize list of features
with ropen(args.list) as fd:
	fts = set( line.rstrip() for line in fd )


#Loop over GFF records
for fields in read_tsv(args.gff, '#'):
	
	#Convert attributes to dictionary
	attrs = attr_to_dict(fields[8])
	
	
	#If this is a child feature of a feature in the list, add it to the list
	if "Parent" in attrs and attrs["Parent"] in fts:
		fts.add( attrs["ID"] )
	
	
	#Report or not report the feature, depending on its status
	if ( args.keep and attrs["ID"] in fts ) or ( not args.keep and attrs["ID"] not in fts ):
		print('\t'.join(fields))
