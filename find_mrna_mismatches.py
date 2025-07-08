#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from os.path import isfile
from Bio.Align import PairwiseAligner

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

parser.add_argument("genome", metavar="GENOME",
	help = "FASTA file referenced by GFF")

parser.add_argument("xcripts", metavar="XCRIPTS",
	help = "FASTA file containing transcript sequences")


#Options
parser.add_argument("-U", "--unstranded",
	help = "Examine both strands",
	action = "store_true")

parser.add_argument("-c", "--count",
	help = "Report number of mismatches, instead of listing mismatches",
	action = "store_true")

parser.add_argument("-A", "--align",
	help = "Report alignments in addition listing occurences",
	action = "store_true")


args = parser.parse_args()



###############################################################################
#Check that all relevant input files exist
require_files( args.gff, args.genome, args.xcripts )

#Will hold number of mismatches found
mismatch_cnt: int = 0


if args.align:
	aln = PairwiseAligner()
	aln.mode = "global"


#Load genes, automatically creating introns with variants
genome = deserialize_fasta(args.genome)
genes = deserialize_gff(args.gff)
xcripts = deserialize_fasta(args.xcripts)
xcripts = truncate_seqids(xcripts)


for gene in genes.values():
	gene.fixup_exons()

for gene in genes.values():
	gene.add_seqs(genome, "gte")


for gene_id in genes.keys():
	#Get sequence of transcripts, based on GFF
	gff_mrna_seq: str = genes[gene_id].transcript.seq
	#Get sequence transcript, based on FASTA
	fa_mrna_seq: str = xcripts[gene_id[:gene_id.index(".path")]]
	
	#Check for matches
	match: bool = gff_mrna_seq in fa_mrna_seq
	
	#If needed, also consider the opposite strand
	if args.unstranded and not match:
		revcomp = reverse_complement(gff_mrna_seq)
		match = revcomp in fa_mrna_seq
	
	
	if not match:
		mismatch_cnt += 1
		if not args.count:
			print(gene_id)
			if args.align:
				if args.unstranded and aln.score(gff_mrna_seq, fa_mrna_seq) < aln.score(revcomp, fa_mrna_seq):
					gff_mrna_seq = revcomp
				
				print(next(aln.align(gff_mrna_seq, fa_mrna_seq)))


#Report number of mismatches if needed
if args.count:
	print(f"{mismatch_cnt} mismatches in {len(genes)} genes")
