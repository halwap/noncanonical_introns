#!/usr/bin/env python3

from introns import *
from argparse import ArgumentParser
from os.path import isfile

################################################################################
#	ARGUMENT PARSING
################################################################################

parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("ANNOT",
	help = "GFF or GTF file contatining annotation")

parser.add_argument("FASTA",
	help = "FASTA file referred to by ANNOT")

parser.add_argument("OUTFILE",
	help = "File to write rectified annotation to")


#Optional arguments: trained models to use
parser.add_argument("-n", "--nonconv",
	help = "Nonconventional model to use instead of default")

parser.add_argument("-c", "--conv",
	help = "Conventional model to use instead of default")


args = parser.parse_args()



################################################################################
#	BODY
################################################################################

if __name__ == "__main__":
	

	#Check that all relevant files exist
	for file in [ args.ANNOT, args.FASTA, args.nonconv, args.conv ]:
		if file is not None and not isfile(file):
			print(f"File {file} does not exist")
			exit()
	
	
	#Reload models, if a specific model was passed
	if args.nonconv:
		models['N'] = load_model(args.nonconv)
	if args.conv:
		models['C'] = load_model(args.conv)
	
	
	#Load genome and genes, creating introns with variants
	genome, genes = create(args.FASTA, args.ANNOT, "gmap")
	
	
	#Get score for all introns of each gene
	predict_all_introns(genes)
	
	
	#Serialize all genes
	#Open file for writing (append mode)
	outfile = open(args.OUTFILE, 'a')
	
	for gene in genes.values():
		gene.finalize_serialize(outfile)
	
	outfile.close()
	

