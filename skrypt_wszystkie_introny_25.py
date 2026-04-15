'''
This file writes all introns from genomes to a file
'''

import os

import src.introns as introns
from collections import defaultdict


path="/mnt/archive/euglena_genomy/"

genom_EG=path+"gracilis/gracilis_dbg2olc.fasta"
geny_EG=path+"gracilis/gracilis_stringtie_strand_informed.gtf"

genom_EH=path+'hiemalis/hiemalis_rascaf.fasta'
geny_EH=path+'hiemalis/stringtie_strand_informed3.gtf'

genom_EL=path+'longa/longa_rascaf.fasta'
geny_EL=path+'longa/longa_stringtie_strand_informed.gtf'

genom_bugtest = '/home/semik/projekty/noncanonical_introns/bugtest.fasta'
geny_bugtest = '/home/semik/projekty/noncanonical_introns/bugtest.gtf'

dlugosc_poczatku_sekwencji=25

for name, genome_eug, genes_eug in [("E.gracilis", genom_EG, geny_EG),
                                    ("E.hiemalis", genom_EH, geny_EH),
                                    ("E.longa", genom_EL, geny_EL)]:#, ("bugtest", genom_bugtest, geny_bugtest)]:
    print(name)
    genome,genes=introns.create_genes(genome_eug, genes_eug)
    print("creating done")

    with open(f'wszystkie_introny_wersje_{name}_{dlugosc_poczatku_sekwencji}zmarginesami.fasta' 'w') as f:
        f.write(";nazwa - nazwagenu;scaffold start - scaffold end;wersja(podst=0);klasa konw;klasa niekonw\n")
        f.write(f";seq - exon[-5:] intron[:{dlugosc_poczatku_sekwencji}] intron[-{dlugosc_poczatku_sekwencji}:] exon[:5]\n")
        print("headline written")
        for nazwa, gene in list(genes.items()):
            #print(nazwa)
            for intron in gene.introns:
                s=intron.sequence
                vconv = intron.is_conventional if intron.is_conventional else 0
                vnonconv = intron.is_nonconventional
                #vstruct = 1 if intron.is_struct_nonconv_1 else 2 if intron.is_struct_nonconv_2 else 0
                f.write(f">{nazwa};{intron.scaffold_start}-{intron.scaffold_end};0;{vconv};{vnonconv}\n")
                f.write(intron.prev_exon.sequence[-5:]+s[:dlugosc_poczatku_sekwencji]+s[-dlugosc_poczatku_sekwencji:]+intron.next_exon.sequence[:5]+"\n")

                for index, var in enumerate(intron.variations):
                    s=var.sequence
                    vconv = var.is_conventional if intron.is_conventional else 0
                    vnonconv = var.is_nonconventional
                    f.write(f">{nazwa};{var.scaffold_start}-{var.scaffold_end};{index};{vconv};{vnonconv}\n")
                    f.write(var.prev_exon.sequence[-5:]+s[:dlugosc_poczatku_sekwencji]+s[-dlugosc_poczatku_sekwencji:]+var.next_exon.sequence[:5]+"\n")

