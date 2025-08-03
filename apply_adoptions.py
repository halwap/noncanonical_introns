#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("pairs", metavar="PAIRS",
					help = "List of pairs of features",
					type = str)
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

#Flags
parser.add_argument("-t", "--type", metavar="TYPE",
					help = "Optional new type for deserted features",
					type = str, default = None)

parser.description = """
Apply adoptions described in PAIRS to GFF.
"""

parser.epilog = """
PAIRS is a 2-column TSV. Field 1 is a feature's ID, field 2 is the ID of
another feature, to which the original feature's children should be moved.
This can be the output of prepare_adoptions.py.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.pairs, args.gff )


if args.type == "":
	raise ValueError(f"Cannot use empty string as new type, use '.' to signify typeless features")


#Child rerouting dictionary
new_parent: dict[str,str] = dict( parse_tsv( args.pairs ) )


for entry in parse_gff(args.gff):

	#If this is a feature to desert, update its type, if requested
	if entry.attrs["ID"] in new_parent and args.type:
		entry.type_ = args.type
	
	
	#If this is a feature to adopt/orphan
	if "Parent" in entry.attrs and entry.attrs["Parent"] in new_parent:
		
		#If this a feature to orphan
		if new_parent[entry.attrs["Parent"]] == '.':
			entry.orphan()
		
		#If this a feature to adopt
		else:
			entry.adopt(new_parent[entry.attrs["Parent"]])
	

	print(entry)