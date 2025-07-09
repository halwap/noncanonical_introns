#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from sequtils import max_orf_len

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file to process")

parser.add_argument("fasta", metavar="FASTA",
	help = "FASTA file referenced by GFF")


#Options
parser.add_argument("-U", "--unstranded",
	help = "Examine both strands",
	default = False, action = "store_true")


args = parser.parse_args()



###############################################################################
#Check that all relevant input files exist
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
