#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from itertools import combinations
from extras import *

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("fasta", metavar="FASTA",
	help = "FASTA file")


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
require_files( args.fasta )


phase("Deserialize")
xcripts = deserialize_fasta(args.fasta)
xcripts = truncate_seqids(xcripts)


phase("Find matching transcripts")
#Get a unique pair of genes in this group
for seqid1, seqid2 in combinations(xcripts, 2):
	#Get the transcript sequences
	t1: str = xcripts[seqid1]
	t2: str = xcripts[seqid2]
	
	#Test if g1 & g2 come from the same transcript
	altseq: bool = are_altseq(t1, t2, args.shared_seq_len)
	
	#If the first test didn't pass, and unstranded analysis was requested, check other
	#orientation
	if not altseq and args.unstranded:
		altseq = are_altseq(reverse_complement(t1), t2, args.shared_seq_len)
	
	#Report the hit if there was one
	if altseq:
		print(f"{seqid1}\t{seqid2}")


phase()