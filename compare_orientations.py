#!/usr/bin/env python3

from argparse import ArgumentParser
from libintrons import *
from os.path import isfile

################################################################################
#   ARGUMENT PARSING
################################################################################

parser = ArgumentParser()


#Positional arguments: input and output files
parser.add_argument("gff", metavar="GFF",
    help = "GFF file contatining annotation")

parser.add_argument("fasta", metavar="FASTA",
    help = "FASTA file referenced by GFF")

parser.add_argument("outfile", metavar="OUTFILE",
    help = "File to write results to")



#Options
parser.add_argument("-n", "--nonconv", metavar="MODEL", required=True,
    help = "Model for scoring nonconventionality")

parser.add_argument("-c", "--conv", metavar="MODEL", required=True,
    help = "Model for scoring conventionality")

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



################################################################################
#   BODY
################################################################################

if __name__ == "__main__":
    

    #Check that all relevant input files exist
    for file in [ args.gff, args.fasta, args.nonconv, args.conv ]:
        if file and not isfile(file):
            raise FileNotFoundError(f"File {file} does not exist")

    
    #If --force was not passed, exit early if output files exist
    if not args.force:
        for file in [ args.outfile ]:
            if file and isfile(file):
                raise FileExistsError(f"File {file} exists, use --force to overwrite")
    
    
    #Load genes, automatically creating introns with variants
    phase("Deserialize")
    nonconv_model = load_model(args.nonconv)
    conv_model = load_model(args.conv)
    genome = deserialize_fasta(args.fasta)
    genes = deserialize_gff(args.gff)
    genes_r = deserialize_gff(args.gff, True)
    
    
    phase("Link exons")
    for gene in genes.values():
        gene.fixup_exons()
    for gene in genes_r.values():
        gene.fixup_exons()
    
    phase("Add introns")
    for gene in genes.values():
        gene.add_introns()
    for gene in genes_r.values():
        gene.add_introns()
        
    phase("Add seqs")
    for gene in genes.values():
        gene.add_seqs(genome)
    for gene in genes_r.values():
        gene.add_seqs(genome)
    
    phase("Add variants")
    for gene in genes.values():
        gene.add_intron_variants(args.min_exon_len)
    for gene in genes_r.values():
        gene.add_intron_variants(args.min_exon_len)
    
    
    #Get score for all introns (& variants) of each gene
    #score_all_introns() automatically triggers new phases
    score_introns(genes.values(), nonconv_model, conv_model,
                  unif_score_from_ss = args.force_conv_variants,
                  batch_size         = args.batch_size,
                  weighted           = not args.unweighted_pairing_scores,
                  pair_scores        = Intron.PairScoresReport(args.pairing_scores_in_scoring))
    score_introns(genes_r.values(), nonconv_model, conv_model,
                  unif_score_from_ss = args.force_conv_variants,
                  batch_size         = args.batch_size,
                  weighted           = not args.unweighted_pairing_scores,
                  pair_scores        = Intron.PairScoresReport(args.pairing_scores_in_scoring))
    
    
    phase("Rectify introns")
    for gene in genes.values():
        gene.rectify_introns()
    for gene in genes_r.values():
        gene.rectify_introns()
    

    phase("Compare orientations")
    with open(args.outfile, 'w') as fd:
        for gene_id in genes.keys():
            gene = genes[gene_id]
            gene_r = genes_r[gene_id]
            if gene.introns:
                avg_score = str( mean( map( attrgetter("unif_score" ), gene.introns ) ))
                avg_score_r = str( mean( map( attrgetter("unif_score" ), gene_r.introns ) ))
            else:
                avg_score = ""
                avg_score_r = ""
            fd.write('\t'.join([gene_id, avg_score, avg_score_r]) + '\n')
    
    
    #End last phase to report its timing
    phase()
