#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "File containing features to keep",
					type = str)

#Options
parser.add_argument("-v", "--invert",
					help = "LIST contains features to discard, rather than keep",
					action = "store_true", default = False)
parser.add_argument("-D", "--deadbeat",
					help = "Children features are kept/discarded independently of their parents",
					action = "store_true", default = False)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )


#Deserialize LIST
fts: set[str] = set( parse_list(args.list) )


#Loop over GFF records
for entry in parse_gff(args.gff):
	
	#Get the ID & Parent attributes of this entry
	id_ = entry.attrs["ID"]
	parent = entry.attrs.get("Parent", "")
	
	
	#If discarded features should have their children orphaned
	if args.deadbeat:
		#If this is a child of a to-be-discarded feature, orphan it
		if ( args.invert and parent in fts ) or ( not args.invert and parent not in fts ):
			entry.orphan()
	
	
	#If discarded features should also have their children discarded
	else:
		#If this is a child of a listed feature list, add it to `fts', so that it meets the same
		#fate as its parent
		if parent != "" and parent in fts:
			fts.add( id_ )
	
	
	#Depending on args.invert and whether this feature is in LIST, keep it or discard it
	if ( args.invert and id_ not in fts ) or ( not args.invert and id_ in fts ):
		print(entry)
