#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from pathlib import Path

###############################################################################

parser = ArgumentParser()

#Positional arguments
parser.add_argument("gff", metavar="GFF",
					help = "GFF file to process",
					type = str)
parser.add_argument("fasta", metavar="FASTA",
					help = "FASTA file referenced by GFF",
					type = str)
parser.add_argument("outfile", metavar="OUTFILE",
					help = "File to write rectified annotation to",
					type = str)
parser.add_argument("stats", metavar="STATS", nargs='?',
					help = "File to write intron statistics to",
					type = str)

#Options
parser.add_argument("-n", "--nonconv", metavar="MODEL",
					help = "Model for scoring nonconventionality",
					type = str, default = Path(__file__).parent / "Models/29_11_K_model.sav")
parser.add_argument("-c", "--conv", metavar="MODEL",
					help = "Model for scoring conventionality",
					type = str, default = Path(__file__).parent / "Models/29_11_NK_model.sav")
parser.add_argument("-f", "--force",
					help = "Force overwriting of generated file(s) if they exist",
					action = "store_true", default = False)
parser.add_argument("-E", "--min-exon-len", metavar="LEN",
					help = "The shortest an exon can become as a result of intron shifting",
					type = int, default = 1)
parser.add_argument("-b", "--batch-size", metavar="SIZE",
					help = "Score introns in batches of this size to reduce peak memory usage",
					type = int, default = 0)
parser.add_argument("-C", "--force-conv-variants",
					help = "Unconditionally prefer variants with conventional splice sites",
					action = "store_true", default = False)
parser.add_argument("-U", "--unweighted-pairing-scores",
					help = "Use unweighted pairing scores in intron scoring",
					action = "store_true", default = False)

args = parser.parse_args()

###############################################################################

#Check that all relevant input files exist
require_files( args.gff, args.fasta, args.nonconv, args.conv )

#If --force was not passed, exit early if output files exist
if not args.force:
	refuse_files( args.stats, args.outfile )


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
#score_introns() automatically triggers new phases
score_introns(genes, nonconv_model, conv_model,
			  unif_score_from_ss = args.force_conv_variants,
			  batch_size         = args.batch_size,
			  weighted           = not args.unweighted_pairing_scores)


phase("Rectify introns")
for gene in genes:
	gene.rectify_introns()


#Re-linking exons and renewing sequences is not necessary at this point in time; they are included
#here, as this script is partly meant as a reference example for using libintrons
#Additionally, calling add_seqs() with "ei" is valuable, as it will verify that the intron
#verification did not alter the sequence of the gene or transcript
#For these reasons, these operations are only carried out if asserts are enabled
if __debug__:
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
		for entry in gene.to_gff():
			fd.write(str(entry))
			fd.write('\n')


if args.stats:
	phase("Intron statistics")
	with wopen(args.stats) as fd:
		#TSV header
		fd.write(to_tsv(*Intron.STATS))
		fd.write('\n')
		
		#TSV body
		for gene in genes:
			for n, intron in enumerate(gene.introns, 1):
				all_stats = intron.get_stats(n)
				for stats in all_stats:
					fd.write(to_tsv(*stats))
					fd.write('\n')


#Call required to report duration of last phase
phase()