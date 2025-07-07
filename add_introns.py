#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from extras import *

################################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

parser.add_argument("fasta", metavar="FASTA",
	help = "FASTA file referenced by GFF")

parser.add_argument("outfile", metavar="OUTFILE",
	help = "File to write annotation with introns to")



#Options
parser.add_argument("-f", "--force",
	help = "Force overwriting of generated file(s) if they exist",
	action = "store_true")


args = parser.parse_args()



################################################################################
#Check that all relevant input files exist
require_files( args.gff, args.fasta )

#If --force was not passed, exit early if output files exist
if not args.force:
	refuse_files( args.outfile )


#Load genes, automatically creating introns with variants
phase("Deserialize")
genome = deserialize_fasta(args.fasta)
genes = deserialize_gff(args.gff)
genes = list(genes.values())


phase("Link exons")
for gene in genes:
	gene.fixup_exons()

phase("Add introns")
for gene in genes:
	gene.add_introns()
	
phase("Add seqs")
for gene in genes:
	gene.add_seqs(genome)

phase("Compute traits")
for gene in genes:
	for intron in gene.introns:
		intron.add_traits()

#Serialize all genes
#Open file for writing
phase("Serialize")
with wopen(args.outfile) as fd:
	for gene in genes:
		gene.serialize(fd)

#End last phase to report its timing
phase()
