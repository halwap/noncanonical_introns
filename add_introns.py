#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *

################################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

parser.add_argument("fasta", metavar="FASTA", nargs='?',
	help = "FASTA file referenced by GFF")


args = parser.parse_args()



################################################################################
#Check that all relevant input files exist
require_files( args.gff, args.fasta )


#Load data
genes = deserialize_gff(args.gff)
genes = list(genes.values())
if args.fasta:
	genome = deserialize_fasta(args.fasta)


#Create introns
for gene in genes:
	gene.fixup_exons()
	#Get introns
	gene.add_introns()

	#If FASTA file was passed, get sequences & splice sites
	if args.fasta:
		gene.add_seqs(genome)
		#Get intron traits, including splice sites
		for intron in gene.introns:
			intron.add_traits()


#Serialize all genes
for gene in genes:
	for entry in gene.to_gff():
		print(entry)
