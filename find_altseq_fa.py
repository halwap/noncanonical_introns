#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *
from itertools import combinations
from extras import are_altseq

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("fasta", metavar="FASTA",
					help = "FASTA file",
					type = str)

#Flags
parser.add_argument("-U", "--unstranded",
					help = "Examine both strands",
					action = "store_true", default = False)
parser.add_argument("-l", "--shared-seq-len", metavar="LEN",
					help = "Minimal shared sequence length to count as altseq hit",
					type = int, default = 1)

parser.description = """
Find pairs of sequences in FASTA which may originate from alternative
sequencing of what is actually the same sequence.
"""

parser.epilog = """
Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.fasta )


xcripts = deserialize_fasta(args.fasta, trunc=True)


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
		print(to_tsv(seqid1,seqid2))
