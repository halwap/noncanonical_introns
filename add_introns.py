#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("fasta", metavar="FASTA", nargs='?',
					help = "FASTA file referenced by GFF",
					type = str)

#Options
parser.add_argument("-E", "--min-exon-len", metavar="LEN",
					help = "The shortest an exon can become as a result of intron shifting",
					type = int, default = 1)

args = parser.parse_args()

###############################################################################

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
		#Get splice site & number of variants
		for intron in gene.introns:
			intron.add_splice_site()
			intron.add_variants(args.min_exon_len)


#Serialize all genes
for gene in genes:
	for entry in gene.to_gff():
		print(entry)
