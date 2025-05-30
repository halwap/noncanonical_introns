from introns import Intron
from collections import defaultdict


def process_file(file_path):
    try:
        with open(file_path) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                yield [int(x) if x.isnumeric() else x for x in line.split()]
    except (IOError, OSError):
        print("Error opening / processing file")


def intron_from_line(line):
    try:
        if isinstance(line[3], str) and line[3] in '+-':
            return Intron(line[0], line[1], line[2], line[3])
        else:
            return Intron(line[0], line[1], line[2])
    except IndexError:
        return Intron(line[0], line[1], line[2])


def extract_introns_from_gtf(file, file_out):
    introns_p = []
    unique = set()
    gene, new_gene, start, end = None, False, 0, 0
    for line in process_file(file):
        if line[2] == 'transcript':
            new_gene = True
            gene = line[9]
        elif line[2] == 'exon':
            if new_gene:
                new_gene = False
                start = line[4]
            else:
                end = line[3] - 1
                scaffold = line[0]
                sign = line[6]
                if not gene:
                    raise Exception()
                i = Intron(scaffold, start, end, gene=gene, strand=sign)
                if ' '.join([scaffold, str(start), str(end)]) not in unique:
                    unique.add(' '.join([scaffold, str(start), str(end)]))
                    introns_p.append(i)
                start = line[4]

    with open(file_out, 'w') as f_out:
        for intron in introns_p:
            f_out.write('\t'.join([str(x) for x in [intron.scaffold, intron.start, intron.end, intron.gene, intron.strand]]))
            f_out.write('\n')


def intron_dict(file):
    my_introns = defaultdict(list)
    for line in process_file(file):
        intron = Intron(line[0], line[1], line[2], gene=line[3], strand=line[4])
        my_introns[intron.scaffold].append(intron)
    return my_introns


def compare_introns(file, introns2, file_out):
    """Compare two sets of introns, classify the way they intersect and save to file the ones that are identical
    in both"""
    all_my_introns = 0
    intron_pairs = []
    no_pair = []
    diff_lengths = []
    for line in process_file(file):
        all_my_introns += 1
        intron = intron_from_line(line)
        pairing = False
        try:
            for i in introns2[intron.scaffold]:
                if intron.intersect(i):
                    intron.gene = i.gene
                    intron_pairs.append([intron, i])
                    pairing = True

            if not pairing:
                no_pair.append(intron)
        except KeyError:
            no_pair.append(intron)

    print('All my introns: ', all_my_introns)
    print('n of pairs: ', len(intron_pairs))
    print('no match: ', len(no_pair))

    exact_match = []
    one_side_match = []
    with open(file_out, 'w') as f_out:
        for i1, i2 in intron_pairs:
            intersection = i1.classify_intersection(i2)
            if intersection == 2:
                exact_match.append([i1, i2])
                f_out.write('\t'.join([i1.scaffold, str(i1.start), str(i1.end)]))
                f_out.write('\n')
            elif intersection == 1:
                one_side_match.append([i1, i2])
                diff_lengths.append(abs(i1.length() - i2.length()))
    return exact_match, one_side_match, diff_lengths
