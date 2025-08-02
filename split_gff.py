#!/usr/bin/env python3

from argparse import ArgumentParser
from formats import *
from extras import bifurcate_gff

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("list", metavar="LIST",
					help = "List of features",
					type = str)
parser.add_argument("out1", metavar="OUTFILE1",
					help = "Outfile for listed features",
					type = str)
parser.add_argument("out2", metavar="OUTFILE2",
					help = "Outfile for listed features",
					type = str)

#Options
parser.add_argument("-f", "--force",
					help = "Force overwriting of generated file(s) if they exist",
					action = "store_true", default = False)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.list )
if not args.force:
	refuse_files( args.out1, args.out2 )


#Deserialize input files
ids: list[str] = list(parse_list( args.list ))
gffs: list[GFF] = list(parse_gff( args.gff ))

#Split feature IDs
ids1, ids2 = bifurcate_gff(gffs, ids)


with wopen(args.out1) as fd1, wopen(args.out2) as fd2:
	for entry in gffs:
		if entry.attrs["ID"] in ids1:
			print(entry, file=fd1)
		if entry.attrs["ID"] in ids2:
			print(entry, file=fd2)
