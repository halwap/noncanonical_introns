#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from libintrons import *
from extras import max_orf_len

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("fasta", metavar="FASTA",
					help = "FASTA file referenced by GFF",
					type = str)

#Flags
parser.add_argument("-U", "--unstranded",
					help = "Examine both strands",
					action = "store_true", default = False)

parser.description = """
Get length of the longest ORF in every 'mRNA' feature in GFF.
"""

parser.epilog = """
GFF is expected to be similar to GMAP's output: 1 'mRNA' feature per 'gene',
>=1 'exon' per 'mRNA'. Other feature types ignored.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.fasta )


genome = deserialize_fasta(args.fasta)
genes = deserialize_gff(args.gff)

#Get transcript sequences
for gene in genes.values():
	gene.fixup_exons()
	gene.add_seqs(genome, "te")


for gene in genes.values():
	#Get transcript sequence
	mrna: str = gene.transcript.seq
	
	
	#Examine ORFs & compare length of longest one
	max_len: int = max_orf_len(mrna)
	
	#Examine reverse ORFs if requested
	if args.unstranded:
		max_len = max( max_len, max_orf_len(reverse_complement(mrna)) )
	
	
	#Report the ORF length
	print(to_tsv(gene.name,max_len))
