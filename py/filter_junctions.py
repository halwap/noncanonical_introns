from intron_comparison import process_file
from introns import Intron
from statistics import median
from itertools import chain
from collections import defaultdict


def _line_to_intron(line):
    """RegTools junctions to Intron."""
    scaffold = line[0]
    start = line[1] + int(line[-2].split(',')[0])
    end = line[2] - int(line[-2].split(',')[1])
    score = line[4]
    return Intron(scaffold, start, end, support=score)


def junctions_to_introns(file_in, file_out):
    """
    Take junctions created by RegTools and save to .bed file only the location and support of introns. Also sort
    the introns by scaffold and then by start for easier analyses later.

    :param file_in: (str) Path to the RegTools file.
    :param file_out: (str) Path to the out file.
    """
    introns_dict = defaultdict(list)
    for line in process_file(file_in):
        i = _line_to_intron(line)
        introns_dict[i.scaffold].append(i)

    with open(file_out, 'w') as f_out:
        for scaff in introns_dict.keys():
            # introns in every scaffold are sorted by start for easier access
            int_list = sorted(introns_dict[scaff], key=lambda intron: intron.start)
            for intron in int_list:
                f_out.write(str(intron))
                f_out.write('\n')


def intron_stats(file):
    """Calculate basic statistics of introns from a file."""
    introns = []
    for line in process_file(file):
        i = Intron(line[0], line[1], line[2], support=line[3])
        introns.append(i)
    print('Number of introns:0', len(introns))
    print('Mean support: ', sum([i.support for i in introns]) / len(introns))
    print('Median support: ', median([intron.support for intron in introns]))


def junctions_to_stats(program, m, path):
    print('Program: ', program)
    file_in = path + 'junctions_' + program + '_m' + m + '.bed'
    file_out = path + 'good_junctions_' + program + '.bed'
    junctions_to_introns(file_in, file_out)
    intron_stats(file_out)


def choose_best_introns(file_in, file_out, cutoff):
    """
    Choose one best intron over every position.

    :param file_in: (str) Path to the .bed file with introns in format: scaffold start end support. All introns from
    a scaffold must come one after another in the file, and within one scaffold introns have to be sorted by start.
    :param file_out: (str) Path to the out file with best introns.
    :param cutoff: (int) Minimum support of the best intron.
    :return: Two dictionaries where key is scaffold and value is the list of introns on the scaffold:
    one containing all the introns from the input file and one with the best introns.
    """
    with open(file_out, 'w') as f_out:

        best_introns = defaultdict(list)
        all_introns = defaultdict(list)

        chrom_old = 'scaffold_0'
        start_old = 0
        end_old = 0
        score_old = 0

        def write_junction():
            junction = '\t'.join([str(x) for x in [chrom, start_old, end_old, score_old]])
            f_out.write(junction)
            f_out.write('\n')
            best_introns[chrom].append(i)

        for line in process_file(file_in):
            chrom, start, end, score = line
            i = Intron(chrom, start, end, support=score)
            all_introns[chrom].append(i)
            # only consider introns with high enough support
            if score < cutoff:
                continue

            if chrom == chrom_old:
                if start < end_old:
                    # still in the same intron
                    if score > score_old:
                        # one best intron in each position
                        start_old, end_old, score_old = start, end, score
                else:
                    # in a new intron, so the old one has to be written down
                    if not start_old - end_old == 0:
                        write_junction()
                    start_old, end_old, score_old = start, end, score

            else:
                # new scaffold, so we need to write down the last intron
                write_junction()
                chrom_old, start_old, end_old, score_old = chrom, start, end, score

        # now we need to write the last one
        write_junction()

        return all_introns, best_introns


def compare_best_introns(all_i, best_i):
    mean_supports = []
    for scaffold in best_i.keys():
        for best_intron in best_i[scaffold]:
            supports = []
            for intron in all_i[scaffold]:
                if best_intron.intersect(intron):
                    supports.append(intron.support)
            mean_supports.append(best_intron.support / sum(supports))
    print('Mean share of the best intron: ', sum(mean_supports) / len(mean_supports))

    sup_of_best = [intron.support for intron in chain.from_iterable(best_i.values())]
    sup_of_all = [intron.support for intron in chain.from_iterable(all_i.values())]
    print('Share of best introns in all: ', sum(sup_of_best) / sum(sup_of_all))
    print('Mean support of the best intron: ', sum(sup_of_all) / len(sup_of_best))
    print('\n')


def introns_for_seq(file_in, file_out, margin):
    """Create a file with introns elongated by a margin to extract sequences with parts of neighbouring exons."""
    with open(file_out, 'w') as f_out:
        for line in process_file(file_in):
            new_line = '\t'.join([str(x) for x in [line[0], line[1] - margin, line[2] + margin]])
            f_out.write(new_line)
            f_out.write('\n')


def introns_for_seq2(file_in, file_out_start, file_out_end, margin1, margin2):
    """Create two files with extracted fragments of margins from both sides of the itron."""
    with open(file_out_start, 'w') as f_out_start:
        with open(file_out_end, 'w') as f_out_end:
            for line in process_file(file_in):
                new_line1 = '\t'.join([str(x) for x in [line[0], line[1] - margin1, line[1] + margin2]])
                f_out_start.write(new_line1)
                f_out_start.write('\n')

                new_line2 = '\t'.join([str(x) for x in [line[0], line[2] - margin2, line[2] + margin1]])
                f_out_end.write(new_line2)
                f_out_end.write('\n')
