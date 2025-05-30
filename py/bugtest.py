import introns
from Bio import SeqIO


def prepare_transcripts(gtf_path, genome_path, fasta_out_path):
    genome = introns.read_genome(genome_path)
    genes = introns.read_genes(gtf_path)
    for_print = []
    for name, gene in genes.items():
        gene.extract_sequence(genome)
        for_print.append([name, gene.transcript.sequence])
    with open(fasta_out_path, 'w') as out_file:
        for gene in for_print:
            out_file.write(">" + gene[0] + '\n')
            out_file.write(gene[1] + '\n')
    return genes


# def get_introns(genes, introns_fasta_out_path):
#     for_print = []
#     for name, gene in genes.items():
#         gene.create_introns()
#         for intron in gene.introns:
#             for_print.append([name, intron.sequence])
#     with open(introns_fasta_out_path, 'w') as out_file:
#         for intron in for_print:
#             out_file.write(">" + intron[0] + '\n')
#             out_file.write(intron[1] + '\n')


if __name__ == '__main__':
    gtf_path = "/Users/pawelhalakuc/robota/genomy/licencjat_julii/noncanonical_introns/bugtest.gtf"
    genome_path = "/Users/pawelhalakuc/robota/genomy/licencjat_julii/noncanonical_introns/bugtest.fasta"
    genes = prepare_transcripts(gtf_path, genome_path, fasta_out_path='bugtest_stringtie.fasta')
    # get_introns(genes, introns_fasta_out_path='bugtest_introns.fasta')

