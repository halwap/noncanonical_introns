#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *
from libintrons import CONV_SS

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)

parser.description = """
Mark each 'intron' feature in GFF with 'splice_site' attribute as conventional
or nonconventional.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

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
