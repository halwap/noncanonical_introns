#!/usr/bin/env python3

from argparse import ArgumentParser
from itertools import combinations
from libintrons import *
from sequtils import are_altseq

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file to process")

parser.add_argument("groups", metavar="GROUPS",
	help = "File with gene groups")

parser.add_argument("fasta", metavar="FASTA",
	help = "FASTA file referenced by GFF")


#Options
parser.add_argument("-U", "--unstranded",
	help = "Examine both strands",
	default = False, action = "store_true")

parser.add_argument("-l", "--shared-seq-len", metavar="LEN",
	help = "Minimal shared sequence length to count as altseq hit",
	default = 1, type = int)


args = parser.parse_args()



###############################################################################
#Check that all relevant input files exist
require_files( args.gff, args.fasta, args.groups )


genome = deserialize_fasta(args.fasta)
genes = deserialize_gff(args.gff)

#Get transcript sequences
for gene in genes.values():
	gene.fixup_exons()
	gene.add_seqs(genome, "te")

#Deserialize gene groups
groups: list[set[str]] = [ set( gene_ids ) for gene_ids in parse_tsv(args.groups) ]


for group in groups:
	#Get a unique pair of genes in this group
	for g1, g2 in combinations(group, 2):
		#Get the transcript sequences
		t1: str = genes[g1].transcript.seq
		t2: str = genes[g2].transcript.seq
		
		#Test if g1 & g2 come from the same transcript
		altseq: bool = are_altseq(t1, t2, args.shared_seq_len)
		
		#If the first test didn't pass, and unstranded analysis was requested, check other
		#orientation
		if not altseq and args.unstranded:
			altseq = are_altseq(reverse_complement(t1), t2, args.shared_seq_len)
		
		#Report the hit if there was one
		if altseq:
			print(to_tsv(g1,g2))


