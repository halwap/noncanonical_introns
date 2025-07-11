#!/usr/bin/env python3

from argparse import ArgumentParser
from itertools import combinations
from libintrons import *
from extras import are_altseq

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("fasta", metavar="FASTA",
					help = "FASTA file referenced by GFF",
					type = str)
parser.add_argument("groups", metavar="GROUPS", nargs='?',
					help = "File with gene groups",
					type = str)

#Flags
parser.add_argument("-U", "--unstranded",
					help = "Examine both strands",
					action = "store_true", default = False)
parser.add_argument("-l", "--shared-seq-len", metavar="LEN",
					help = "Minimal shared sequence length to count as altseq hit",
					type = int, default = 1)

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.fasta, args.groups )


genome = deserialize_fasta(args.fasta)
genes = deserialize_gff(args.gff)

#Get transcript sequences
for gene in genes.values():
	gene.fixup_exons()
	gene.add_seqs(genome, "te")


#Deserialize gene groups if it was passed, or put all genes in the GFF into a single group
if args.groups:
	groups: list[list[str]] = list[ parse_tsv(args.groups) ]
else:
	groups: list[Iterable[str]] = [ genes.keys() ]


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


