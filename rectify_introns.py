#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from os.path import realpath, dirname
from extras import *

###############################################################################
parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
	help = "GFF file contatining annotation")

parser.add_argument("fasta", metavar="FASTA",
	help = "FASTA file referenced by GFF")

parser.add_argument("outfile", metavar="OUTFILE",
	help = "File to write rectified annotation to")

parser.add_argument("stats", metavar="STATS", nargs='?',
	help = "File to write intron statistics to")



#Options
parser.add_argument("-n", "--nonconv", metavar="MODEL",
	help = "Model for scoring nonconventionality",
	default = dirname(realpath(__file__)) + "/Models/29_11_K_model.sav")

parser.add_argument("-c", "--conv", metavar="MODEL",
	help = "Model for scoring conventionality",
	default = dirname(realpath(__file__)) + "/Models/29_11_NK_model.sav")

parser.add_argument("-f", "--force",
	help = "Force overwriting of generated file(s) if they exist",
	action = "store_true")

parser.add_argument("-E", "--min-exon-len", metavar="LEN",
	help = "The shortest an exon can become as a result of intron shifting",
	default = 1, type = int)

parser.add_argument("-b", "--batch-size", metavar="SIZE",
	help = "Score introns in batches of this size to reduce peak memory usage",
	default = 0, type = int)

parser.add_argument("-C", "--force-conv-variants",
	help = "Unconditionally prefer variants with conventional splice sites",
	default = False, action = "store_true")

parser.add_argument("-U", "--unweighted-pairing-scores",
	help = "Use unweighted pairing scores in intron scoring",
	default = False, action = "store_true")

parser.add_argument("-P", "--pairing-scores-in-scoring", metavar="HOW",
	help = "How to treat pairing scores during intron scoring",
	default = Intron.PairScoresReport.KEEP_ALL.value, type = int)

parser.epilog = """
Valid values of HOW are:
0 - Take as-is,
1 - Drop all,
2 - Copy highest,
3 - Drop all but highest
"""


args = parser.parse_args()



###############################################################################
#Check that all relevant input files exist
require_files( [ args.gff, args.fasta, args.nonconv, args.conv ] )

#If --force was not passed, exit early if output files exist
if not args.force:
	refuse_files( [ args.stats, args.outfile ] )


#Load genes, automatically creating introns with variants
phase("Deserialize")
nonconv_model = load_model(args.nonconv)
conv_model = load_model(args.conv)
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

phase("Add variants")
for gene in genes:
	gene.add_intron_variants(args.min_exon_len)


#Get score for all introns (& variants) of each gene
#score_all_introns() automatically triggers new phases
score_introns(genes, nonconv_model, conv_model,
			  unif_score_from_ss = args.force_conv_variants,
			  batch_size         = args.batch_size,
			  weighted           = not args.unweighted_pairing_scores,
			  pair_scores        = Intron.PairScoresReport(args.pairing_scores_in_scoring))


phase("Rectify introns")
for gene in genes:
	gene.rectify_introns()

phase("Re-link exons")
for gene in genes:
	gene.fixup_exons()

phase("Re-add seqs")
for gene in genes:
	gene.add_seqs(genome, "ei")


#Serialize all genes
#Open file for writing
phase("Serialize")
with wopen(args.outfile) as fd:
	for gene in genes:
		gene.serialize(fd)


if args.stats:
	phase("Intron statistics")
	with wopen(args.stats) as fd:
		#TSV header
		fd.write(
			'\t'.join([
				"gene",
				"transcript",
				"intron_idx",
				"variant_cnt",
				"variant_rank",
				"scaffold",
				"strand",
				"start",
				"end",
				"unif_score",
				"c_score",
				"nc_score",
				"score_range",
				"score_outpace",
				"splice_site",
				"splice_site_is_conv",
				"e-5",
				"i10",
				"i-10",
				"e5",
				"prev_exon_y",
				"start_r",
				"end_y",
				"next_exon_r",
				"cagctg",
				"pairing_score_1",
				"pairing_score_2",
				"pairing_score_3",
				"pair_3_-6",
				"pair_4_-7",
				"pair_5_-8"
			]) + '\n'
		)
		
		#TSV body
		for gene in genes:
			for n, intron in enumerate(gene.introns, 1):
				intron.write_stats(fd, n)


#End last phase to report its timing
phase()
