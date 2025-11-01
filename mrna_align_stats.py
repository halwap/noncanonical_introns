#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from libintrons import *
from Bio.Align import PairwiseAligner
from extras import get_aln_layout

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("genome", metavar="GENOME",
					help = "FASTA file referenced by GFF",
					type = str)
parser.add_argument("xcripts", metavar="XCRIPTS",
					help = "FASTA file containing transcript sequences",
					type = str)

#Flags
parser.add_argument("-U", "--unstranded",
					help = "Examine both strands",
					action = "store_true", default = False)

parser.description = """
Given GFF output by GMAP as a result of mapping XCRIPTS to FASTA, align each
'mRNA' feature to its original transcript, and return the alignment statistics.
"""

parser.epilog = """
GFF is expected to be GMAP output, or very closely similar to it in format.

Any path can be '-' to read from stdin. Writes to stdout.
"""

args = parser.parse_args()

###############################################################################

require_files( args.gff, args.genome, args.xcripts )


aln = PairwiseAligner()
aln.mode = "global"
#Disable penalty for leading & trailing gaps
aln.left_gap_score=0
aln.right_gap_score=0


#Load genes, automatically creating introns with variants
genome = deserialize_fasta(args.genome)
genes = deserialize_gff(args.gff)
xcripts = deserialize_fasta(args.xcripts, trunc=True)

for gene in genes.values():
	gene.fixup_exons()

for gene in genes.values():
	gene.add_seqs(genome, "gte")


#Print TSV header
print(to_tsv(
		"id",				#GFF feature ID
		"base_seq_len",		#Length of base sequence
		"map_seq_len",		#Length of mapped sequence
		"aln_len",			#Length of alignment
		"matches",
		"mismatches",
		"gaps"
	))


for gene_id in genes.keys():
	#Get sequence of transcripts, based on GFF
	gff_seq: str = genes[gene_id].transcript.seq
	#Get sequence transcript, based on FASTA
	fa_seq: str = xcripts[gene_id[:gene_id.index(".path")]]
	
	
	#Compute alignment
	if args.unstranded:
		#If both strands are to be examined, get the configuration which yields
		#a higher score
		revcomp = reverse_complement(gff_seq)
		if aln.score(gff_seq, fa_seq) >= aln.score(revcomp, fa_seq):
			alignment = next(aln.align(gff_seq, fa_seq))
		else:
			alignment = next(aln.align(revcomp, fa_seq))
	else:
		#Otherwise, just get the alignment as-is
		alignment = next(aln.align(gff_seq, fa_seq))
	
	
	#Get layout of alignment, stripping out terminal gaps
	layout = get_aln_layout(alignment)
	layout = layout.strip()


	#Raport
	print(to_tsv(
		gene_id,
		len(fa_seq),
		len(gff_seq),
		alignment.length,
		layout.count('|'),
		layout.count('.'),
		layout.count(' ')
	))
