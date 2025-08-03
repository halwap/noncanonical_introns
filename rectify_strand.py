#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from formats import *

###############################################################################

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

#Positional arguments
parser.add_argument("fasta", metavar="FASTA",
					help = "FASTA file containing transcript sequences",
					type = str)
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("blastx", metavar="BLASTX",
					help = "Two-column BLASTX results",
					type = str)
parser.add_argument("orfs", metavar="ORFS",
					help = "TD2.LongOrfs GFF output",
					type = str)
parser.add_argument("scores", metavar="SCORES",
					help = "Three-column mapping intron scores",
					type = str)

parser.description = """
Determine the best orientation for every transcript, and apply this orientation
to transcript mappings in GFF.
"""

parser.epilog = """
GFF should be GMAP output, or of similar format to it.

BLASTX is a 2-column TSV with BLASTX results, matching BLAST's
'6 qseqid qstrand' output format specifier.

ORFS is 'longest_orfs.gff3' as output by TD2, filtered to 'CDS' features only.

SCORES is a 3-column TSV. Field 1 is a 'gene' feature ID, fields 2 & 3 are its
mean intron scores when scored as-is and in inverted orientation, respectively.

Any path can be '-' to read from stdin. Writes rectified GFF to stdout, and
rectification statistics to stderr.
"""

args = parser.parse_args()

###############################################################################

require_files( args.fasta, args.gff, args.blastx, args.orfs, args.scores )


#Get xcript IDs, truncated to only the first word
xcript_ids = list( deserialize_fasta(args.fasta, trunc=True).keys() )


#Will hold the decisions for each transcript, based on BLASTX results
blastx_res: dict[str:str] = { xcript_id:'.' for xcript_id in xcript_ids }

#Loop over BLASTX results, to fill out `blastx_red' properly
for xcript_id, strand in parse_tsv(args.blastx):
	blastx_res[xcript_id] = strand


#Will hold the decisions for each transcript, based on ORF lengths
orf_res: dict[str:str] = { xcript_id:'.' for xcript_id in xcript_ids }
#Will hold the length of the longest ORF found for a given transcript
orf_lens: dict[str:int] = { xcript_id:0 for xcript_id in xcript_ids }

#Fill out `orf_res' properly
for entry in parse_gff(args.orfs):
	if len(entry) > orf_lens[entry.seqid]:
		orf_lens[entry.seqid] = len(entry)
		orf_res[entry.seqid] = entry.strand


#Will hold the decisions for each transcript, based on intron scores
#These score are actually per transcript mapping, of which there may be
#>1 per transcript; each mapping is parsed independently, and the per-mapping
#decisions are flattened to be per-transcript, if needed
scores_res: dict[str:list[str]] = { xcript_id:[] for xcript_id in xcript_ids }

#Fill out `scores_res'
for id_, score1, score2 in parse_tsv(args.scores):
	#Get transcript id
	xcript_id = id_[:id_.index('.path')]
	
	score1 = float(score1)
	score2 = float(score2)

	#Current orientation is better
	if score1 > score2:
		scores_res[xcript_id].append('+')
	#Reversed orientation is better
	elif score2 > score1:
		scores_res[xcript_id].append('-')
	#Tie
	else:
		scores_res[xcript_id].append('.')


#Flatten per-mapping decisions to per-transcript
for xcript_id, dec in scores_res.items():
	#If the transcript has no mappings, skip it
	if len(dec) == 0:
		scores_res[xcript_id] = '.'

	#If the transcript has one mapping, take it as the per-transcript decision
	elif len(dec) == 1:
		scores_res[xcript_id] = dec[0]

	#If the transcript has multiple mappings, flatten their decisions to one
	else:
		#Did any mapping get a '+'?
		plus = '+' in dec
		#Did any mapping get a '-'?
		minus = '-' in dec
		
		
		#If some mappings got a '+' and some a '-', skip the transcript
		if plus and minus:
			scores_res[xcript_id] = '.'
		
		#If mappings only ever got a '+' or a '.', take '+'
		elif plus:
			scores_res[xcript_id] = '+'
		
		#Likewise, if mappings only ever got a '-' or a '.', take '-'
		elif minus:
			scores_res[xcript_id] = '-'
		
		#Otherwise, mappings only ever got a '.'; take it
		else:
			scores_res[xcript_id] = '.'


#Flatten each transcript's 3 decisions into one,
#with priority BLASTX > ORF length > intron scores
final_res: dict[str:str] = { xcript_id:'.' for xcript_id in xcript_ids }
for xcript_id in xcript_ids:
	blastx = blastx_res[xcript_id]
	orf    = orf_res[xcript_id]
	scores = scores_res[xcript_id]
	
	#BLASTX results have highest priority
	if blastx != '.':
		final_res[xcript_id] = blastx
	
	#ORF lengths take second priority
	elif orf != '.':
		final_res[xcript_id] = orf
	
	#compare_orientations.py have lowest priority
	elif scores != '.':
		final_res[xcript_id] = scores


#Write TSV listing the decisions for each transcript
eprint(to_tsv( "xcript_id", "blastx", "orf_len", "intron_score", "decision" ))

for xcript_id in xcript_ids:
	eprint(to_tsv(
		xcript_id,
		blastx_res[xcript_id],
		orf_res[xcript_id],
		scores_res[xcript_id],
		final_res[xcript_id]
	))


#Finally, process GFF
for entry in parse_gff(args.gff):
	#Get transcript id
	xcript_id = entry.attrs["Name"]
	assert xcript_id in final_res
	
	#Invert strand if necessary
	if final_res[xcript_id] == '-':
		if entry.strand == '+':
			entry.strand = '-'
		elif entry.strand == '-':
			entry.strand = '+'
	
	#Print feature
	print(entry)


