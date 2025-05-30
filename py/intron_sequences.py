from introns import Intron
import re
import collections, functools, operator
from collections import defaultdict


def file_to_seq_introns(file, margin):
    """Extract introns with sequences from a fasta file"""
    introns = []
    with open(file) as f_in:
        for line in f_in.readlines():
            if line[0] == '>':
                line = re.split('[>:\-]', line.strip())
                scaffold, start, end = line[1:]
            else:
                i = Intron(scaffold, int(start), int(end), margin_left=margin, margin_right=margin, sequence=line.strip())
                introns.append(i)
    return introns


def oligofreq(sequence, k):
    """Calculate the occurence of oligonucleotides of length n"""
    vector = defaultdict(int)
    for pos in range(len(sequence) - k + 1):
        vector[sequence[pos : pos + k]] += 1
    # return np.array([vector["".join(kmer)] for kmer in product("ACGT", repeat=k)])
    return vector


def seq_statistics(introns):
    """Calculate basic statistics of biological sequences"""
    length = sum([intron.length() for intron in introns]) / len(introns)
    print('Number of introns: ', len(introns))
    print('Mean length: ', length)

    def gc_content(intron):
        gc = 0
        for i, nucleotide in enumerate(intron.sequence):
            if nucleotide in 'GC' and i in range(intron.margin_left, intron.length() + intron.margin_left):
                gc += 1
        return gc / (intron.length())

    gc = sum([gc_content(intron) for intron in introns]) / len(introns)
    print('Mean gc: ', gc)

    print('Tetranucleotides: ')
    ini_dict = []
    for intron in introns:
        ini_dict.append(oligofreq(intron.sequence, 4))
    result = dict(functools.reduce(operator.add, map(collections.Counter, ini_dict)))
    print([(j, v) for j, v in sorted(result.items(), key=lambda item: item[1], reverse=True)][:10])
#     return dict([(j, v) for j, v in sorted(result.items(), key=lambda item: item[1], reverse=True)][:10])

    def extract_junction(intron):
        return intron.sequence[intron.margin_left:intron.margin_left + 2] +\
               intron.sequence[-intron.margin_right - 2:-intron.margin_right]

    junctions = defaultdict(int)
    for intron in introns:
        junction = extract_junction(intron)
        junctions[junction] += 1

    top10 = [(j, v) for j, v in sorted(junctions.items(), key=lambda item: item[1], reverse=True)][:10]
    print('Top junctions: ', top10)
    junctions = [(j, v) for j, v in sorted(junctions.items(), key=lambda item: item[1], reverse=True)]
    return junctions
