from collections import defaultdict
from copy import copy
from re import search
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import subprocess

nonconventional_model = pickle.load(open('29_11_K_model.sav', 'rb'))
conventional_model = pickle.load(open('29_11_NK_model.sav', 'rb'))
#nonconventional_model = pickle.load(open('./files_for_classifiers/29_11_NK_model.sav', 'rb'))
#conventional_model = pickle.load(open('./files_for_classifiers/29_11_K_model.sav', 'rb'))
#nonconventional_model = pickle.load(open('./files_for_classifiers/15_11_NK_model.sav', 'rb'))
#conventional_model = pickle.load(open('./files_for_classifiers/15_11_K_model.sav', 'rb'))

def get_conventional_model(): return conventional_model
def get_nonconventional_model(): return nonconventional_model


def getter_setter_gen(name, type_):
    def getter(self):
        return getattr(self, "__" + name)

    def setter(self, value):
        if not isinstance(value, type_) and value is not None:
            raise TypeError("%s attribute must be set to an instance of %s" % (name, type_))
        setattr(self, "__" + name, value)
    return property(getter, setter)


def auto_attr_check(cls):
    new_dct = {}
    for key, value in cls.__dict__.items():
        if isinstance(value, type):
            value = getter_setter_gen(key, value)
        new_dct[key] = value
    # Creates a new class, using the modified dictionary as the class dict:
    return type(cls)(cls.__name__, cls.__bases__, new_dct)


class GenomicSequence:
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None):
        self.scaffold_name = scaffold_name
        self.scaffold_start = scaffold_start
        self.scaffold_end = scaffold_end
        if sequence and len(sequence) != self.length():
            raise ValueError('Incorrect sequence length.')
        self.sequence = sequence
        if strand and strand not in ['+', '-', '.']:
            raise ValueError('Strand can only be +, - or .')
        else:
            self.strand = strand

    def length(self):
        return self.scaffold_end - self.scaffold_start
    
    def __repr__(self):
        return ' '.join([self.scaffold_name, str(self.scaffold_start), str(self.scaffold_end)])
        
    def __str__(self):
        return self.__repr__()


class Gene(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='', transcript=None, exons=None, introns=None,
                 name=''):
        # if strand == '-':
        #     start, end = end, start
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.transcript = transcript
        self.exons = [] if exons is None else exons
        self.introns = [] if introns is None else introns
        self.name = name
        self.expanded_sequence = ''
        self.expansion_left = 0
        self.expansion_right = 0
        self.introns_dict = {}
        self.coverage = coverage
        
    def append_exons(self, exon):
            self.exons.append(exon)
            
    def add_exons(self, genes_data_type):
        if self.strand == '+':
            self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
        elif self.strand == '-':
            self.working_exons.sort(key=lambda _exon: _exon.scaffold_start, reverse=True)
        else:
            raise ValueError('gene strand not in [+, -]'.format(self.strand))
        prev = None
        for exon in self.working_exons:
            if prev:
                if genes_data_type == 'manual':
                    if prev.scaffold_end == exon.scaffold_start:
                        prev.merge_exons(exon)
                        continue
                prev.next_exon = exon
                exon.prev_exon = prev
            self.exons.append(exon)
    
    def append_introns(self, intron):
        self.introns.append(intron)
    
    def create_intron_structures(self, path_working="./ViennaRNA-2.4.18/working_fastas/", path_RNAfold = '/usr/local/bin/'):
        
        p=plik.split('/')[-1].split('.')[0]
        command = [path_RNAfold+'RNAfold', '--noPS',  '--command=constr.txt', '-j8', '-i', path_working+p+'.fasta']
        with open(path_working+p+'.dbn', 'w') as output_file:
            _p = subprocess.run(command, stdout=output_file)

    def add_exons(self, genes_data_type):
        if self.strand == '+':
            self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
        elif self.strand == '-':
            self.working_exons.sort(key=lambda _exon: _exon.scaffold_start, reverse=True)
        else:
            raise ValueError('gene strand not in [+, -]'.format(self.strand))
        prev = None
        for exon in self.working_exons:
            if prev:
                if genes_data_type == 'manual':
                    if prev.scaffold_end == exon.scaffold_start:
                        prev.merge_exons(exon)
                        continue
                prev.next_exon = exon
                exon.prev_exon = prev
            self.exons.append(exon)
            prev = exon

    def extract_sequence(self, genome):
        # elif self.strand == '-':
        #     sequence = genome[self.scaffold_name][self.scaffold_end:self.scaffold_start]
        # else:
        #     raise Exception('cojest')
        scaffold_seq = genome[self.scaffold_name]
        sequence = scaffold_seq[self.scaffold_start:self.scaffold_end]
        expanded_sequence, expansion_left, expansion_right = self.get_expanded_sequence(scaffold_seq)
        if self.strand == '-':
            self.sequence = reverse_complement(sequence)
            self.expanded_sequence = reverse_complement(expanded_sequence)
            self.expansion_left = expansion_right
            self.expansion_right = expansion_left
        else:
            self.sequence = sequence
            self.expanded_sequence = expanded_sequence
            self.expansion_right = expansion_right
            self.expansion_left = expansion_left
        for exon in self.exons:
            if self.strand == '+':
                exon.sequence = self.sequence[exon.scaffold_start - self.scaffold_start:exon.scaffold_end - self.scaffold_start]
            elif self.strand == '-':
                exon.sequence = self.sequence[- exon.scaffold_end + self.scaffold_end:- exon.scaffold_start + self.scaffold_end]
            else:
                print(self.name, exon.scaffold_name, exon.scaffold_start, exon.scaffold_end)
                # raise Exception('co jest')
        transcript_sequence = self.get_transcript_sequence()
        self.transcript = Transcript(self.scaffold_name, self.scaffold_start, self.scaffold_end, strand=self.strand,
                                     sequence=transcript_sequence)

    def get_expanded_sequence(self, scaffold_seq):
        start = max(0, self.scaffold_start - 500)
        end = min(len(scaffold_seq), self.scaffold_end + 500)
        expanded_sequence = scaffold_seq[start:end]
        expansion_left = self.scaffold_start - start
        expansion_right = end - self.scaffold_end
        return expanded_sequence, expansion_left, expansion_right

    def get_transcript_sequence(self):
        exons_seqs = []
        for exon in self.exons:
            exons_seqs.append((exon.sequence, exon.scaffold_start))
        if self.strand == '+':
            exons_seqs.sort(key=lambda tup: tup[1])
        elif self.strand == '-':
            exons_seqs.sort(key=lambda tup: tup[1], reverse=True)
        sequence = ''.join([exon_seq[0] for exon_seq in exons_seqs])
        return sequence

    def get_transcript_with_gaps_sequence(self, expanded=False):
        to_be_joined = []
        start, end = None, None
        reverse = True if self.strand == '-' else False
        exons_sorted = sorted(self.exons, key=lambda obj: obj.start, reverse=reverse)
        if not reverse:
            for exon in exons_sorted:
                start = exon.scaffold_start
                if end:
                    intron_length = start - end
                    assert intron_length > -1
                    to_be_joined.append((''.join(['-' for i in range(intron_length)])))
                to_be_joined.append(exon.sequence)
                end = exon.scaffold_end
        else:
            for exon in exons_sorted:
                start = exon.scaffold_end
                if end:
                    intron_length = end - start
                    assert intron_length > -1
                    to_be_joined.append((''.join(['-' for i in range(intron_length)])))
                to_be_joined.append(exon.sequence)
                end = exon.scaffold_start
        sequence = ''.join(to_be_joined)
        if expanded:
            sequence = self.expansion_left * '-' + sequence + self.expansion_right * '-'
        return sequence

    def create_introns(self):
        if len(self.exons) < 1:
            raise ValueError('No exons specified.')
            return
        else:
            self.introns = []
            start, end = 0, 0
            for exon in self.exons:
                prev_e = exon.prev_exon
                end = exon.scaffold_start
                if start:
                    if self.strand == '+':
                        sequence = self.sequence[start - self.scaffold_start:end - self.scaffold_start]
                    elif self.strand == '-':
                        sequence = self.sequence[- end + self.scaffold_end: - start + self.scaffold_end]
                    intr=Intron(self.scaffold_name, scaffold_start = start, scaffold_end = end, strand=self.strand, sequence=sequence, gene=self, prev_exon=prev_e, next_exon=exon)
                    self.append_introns(intr)
                    
                    prev_e.next_intron = intr
                    exon.prev_intron=intr

                    prev_e=exon
                    # elif self.strand == '-':
                    #     self.append_introns(Intron(self.scaffold_name, end, start))
                    # else:
                    #     raise Exception('co jest')
                start = exon.scaffold_end
        for intron in self.introns:
            intron.movable_boundary_no_margins()
            intron.conventional_version()
            intron.nonconventional_version()
            intron.noncanonical_structural_versions()# za wolno działa albo nie umiem
        
        # for intron in self.introns:
        #     if self.strand == '-':
        #         print(str(intron), intron.scaffold_start, self.scaffold_start, intron.scaffold_end, self.scaffold_start, self.strand, intron.strand)
        #         intron.sequence = self.sequence[intron.scaffold_start - self.scaffold_start:intron.scaffold_end - self.scaffold_start]
        #         print(self.sequence)
        #         print(intron.sequence)
        #     elif self.strand == '+':
        #         print(str(intron), intron.scaffold_end, self.scaffold_end, intron.sscaffold_tart, self.scaffold_end, self.strand, intron.strand)
        #         intron.sequence = self.sequence[- intron.scaffold_end + self.scaffold_end:- intron.scaffold_start + self.scaffold_end]
        #         print(self.sequence)
        #         print(intron.sequence)


@auto_attr_check
class Intron(GenomicSequence):
    """
    This is a class for working with biological introns.

    :param scaffold_name: (str) Name of scaffold or chromosome on which the intron is located.
    :param start: (int) Location of the first base of the intron (within scaffold).
    :param end: (int) Location of the first base of the next exon (within scaffold).
    :param start: (int) Location of the first base of the intron (within gene).
    :param end: (int) Location of the first base of the next exon (within gene).
    :param gene: (str) Gene in which the intron is located.
    :param strand (str): Optional, defines the strand on which intron is located. Must be either + or -.
    :param support: (int) Optional, how many reads support the intron.
    :param margin_left: (int) Optional, how many nucleotides from the preceding exon are included.
    :param margin_right: (int) Optional, how many nucleotides from the following exon are included.
    :param margin_left_seq: (str) Optional, end sequence from the preceding exon.
    :param margin_right_seq: (str) Optional, beginning sequence of the following exon.
    :param sequence: (str) Optional, genomic sequence of the intron.
    :param is_conventional: (int) Optional, number of the conventional class the intron belongs to; if isn't conventional then 0.
    :param is_nonconventional: (int) Optional, number of the nonconventional class the intron belongs to; if isn't nonconventional then 0.
    :param best_conv_var: (int) Optional, unique to main intron, absent in variations; number of the best conventional class out of all variations.
    :param best_nonconv_var: (int) Optional, unique to main intron, absent in variations; number of the best nonconventional class out of all variations.
    """
    scaffold_name = str
    scaffold_start = int
    scaffold_end = int
    gene_end = int
    gene_start = int
    gene = Gene
    strand = str
    support = int
    margin_left = int
    margin_right = int
    sequence = str
    margin_left_seq = str
    margin_right_seq = str
    variations = list
    is_conventional = int
    is_nonconventional = int
    variations_struct_1 = int
    variations_struct_2 = int
    is_struct_nonconv_1 = int
    is_struct_nonconv_2 = int
    man_annotation = str
    test_annotation = str

    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None, gene=None, support=None, margin_left=0,
                 margin_right=0, margin_left_seq='', margin_right_seq='', prev_exon=None, next_exon=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.gene = gene
        self.support = support
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_left_seq = margin_left_seq
        self.margin_right_seq = margin_right_seq
        self.variations = []
        self.is_conventional = 0
        self.is_nonconventional = 0
        self.prev_exon = prev_exon
        self.next_exon = next_exon
        self.best_conv_var = 0  # variation of the intron with the best conventional version
        self.best_nonconv_var = 0  # variation of the intron with the best nonconventional version
        self.best_conv_obj = None
        self.best_nonconv_obj = None
        self.variations_struct_1 = None
        self.variations_struct_2 = None
        self.is_struct_nonconv_1 = None
        self.is_struct_nonconv_2 = None
        self.man_annotation = man_annotation
        self.man_variant = None
        self.test_annotation = test_annotation
        self.test_global_annotation = test_global_annotation
        self.test_score = test_score
        self.test_score_max = -1000
        self.test_score_min = 1000
        self.test_best_k_var = self
        self.test_best_nk_var = self
        self.test_k = False
        self.test_nk = False
        self.canonical_borders = False
        self.gc_content = gc_content
        self.polypyrimidine_tract = polypyrimidine_tract
        self.conserved_pairing_score = 0
        self.ML_characteristic = np.zeros(8)
        self.ML_class = None #[is k, is nk]
        self.ML_is_nonconventional = False
        self.ML_best_conv_version = None
        self.ML_best_nonconv_version = None
        self.ML_conv_score = 0
        self.ML_nonconv_score = 0
        # TODO przy zmienianiu podstawowego intronu trzeba przepisac best_(non)conv_var i liste wariacji
        #      - warianty ich nie maja i zawsze maja nie miec
        if self.strand == '-':
            self.gene_start, self.gene_end = -self.scaffold_end + self.gene.scaffold_end,\
                                             -self.scaffold_start + self.gene.scaffold_end
        elif self.strand == '+':
            self.gene_start, self.gene_end = self.scaffold_start - self.gene.scaffold_start,\
                                             self.scaffold_end - self.gene.scaffold_start

    def movable_boundary_no_margins(self):
        # moglem zepsuc
        # ja zaraz jeszcze bardziej popsuje
        """
        Check if there are repeats on the intron junctions so the intron position could be shifted without changing
        transcript sequence. If there are, add new possible introns to self.variations.
        """
        
        if not self.prev_exon or not self.next_exon:
            return
        mls, mrs = self.prev_exon.sequence, self.next_exon.sequence
        
        i = 1
        # start checking for repeats left from the junction
        check = 'left'
        
        while True:
            new_prev_exon = copy(self.prev_exon)
            new_next_exon = copy(self.next_exon)
            if check == 'left': #checking to the left
                if i > len(mls) or i > len(self.sequence) or mls[-i] != self.sequence[-i]:
                    check = 'right'
                    i = 0
                    continue
                new_seq=mls[-i:]+self.sequence[:-i]
                new_prev_exon.sequence = new_prev_exon.sequence[:-i]
                new_prev_exon.scaffold_end = new_prev_exon.scaffold_end-i
                new_next_exon.sequence = self.sequence[-i:]+new_next_exon.sequence
                new_next_exon.scaffold_start = new_next_exon.scaffold_start-i
                #creating new variation moved to the left
                new_variation = Intron(self.scaffold_name, scaffold_start=self.scaffold_start - i, scaffold_end=self.scaffold_end - i, \
                                       gene=self.gene, sequence=new_seq, prev_exon=new_prev_exon, next_exon=new_next_exon)
                
            else: #checking to the right
                if i + 1 > len(mrs) or i + 1 > len(self.sequence) or self.sequence[i] != mrs[i]:
                    break
                new_seq=self.sequence[i:]+mrs[:i]
                new_prev_exon.sequence = new_prev_exon.sequence+self.sequence[:i]
                new_prev_exon.scaffold_end = new_prev_exon.scaffold_end+i
                new_next_exon.sequence = new_next_exon.sequence[i:]
                new_next_exon.scaffold_start = new_next_exon.scaffold_start+i
                #creating new variation moved to the right
                new_variation = Intron(self.scaffold_name, scaffold_start=self.scaffold_start + i, scaffold_end=self.scaffold_end + i, \
                                       gene=self.gene, sequence=new_seq, prev_exon=new_prev_exon, next_exon=new_next_exon)
            self.variations.append(new_variation)
            i += 1

    def check_conventional(self):
        """ Check if the intron junctions suggest the intron is conventional."""
        left_anchor = self.sequence[0:2]
        right_anchor = self.sequence[-2:]
        
        if left_anchor in ['GT', 'GC'] and right_anchor == 'AG':
            return True
        elif left_anchor == 'CT' and right_anchor in ['AC', 'GC']:
            return True
        else:
            return False

    def check_nonconventional(self):
        """Check if intron may be nonconventional according to our current knowledge, meaning it can form
        secondary structure in specific positions. Also check if variations with shifted junctions may be
        nonconventional."""
        #funkcje complimentary wyciagnelam poza te funkcje

        left_anchor = self.sequence[3: 5]
        right_anchor = self.sequence[-7:-5]

        # try:
        if complimentary(left_anchor[0], right_anchor[1]) and complimentary(left_anchor[1], right_anchor[0]):
            return True
        else:
            # checking variations
            if len(self.variations) > 0:
                for son in self.variations:
                    if son.check_nonconventional():
                        return True
            return False
            
    def conventional_version(self):
        if self.sequence[0:2] in ['GT', 'GC'] and self.sequence[-2:] == 'AG':
            i = 4
        else:
            return
        
        if not self.prev_exon or not self.next_exon:
            return
        mls, mrs = self.prev_exon.sequence[-3:], self.next_exon.sequence[:3]
        
        #mls, mrs = self.margin_left_seq, self.margin_right_seq
        #if check == 0:
        if self.sequence[-3] == 'C':
            i = 3
            if mls and mls[-1] == 'G' and \
                mrs and mrs[0] == 'G':
                i = 2
                if len(mls) > 1 and len(mrs) > 1 and \
                    mls[-2] == 'A' and mrs[1] == 'T':
                    i = 1
        if self.sequence[1] == 'C': i += 4
        self.is_conventional = i
        # else:  # if check = 1
        #     if self.sequence[2] == 'G':
        #         i = 3
        #         if self.margin_left_seq and self.margin_right_seq and \
        #            self.margin_left_seq[-1] == 'C' and self.margin_right_seq[0] == 'C':
        #             i = 2
        #             if len(self.margin_left_seq) > 1 and len(self.margin_right_seq) > 1 and \
        #                self.margin_left_seq[-2] == 'A' and self.margin_right_seq[1] == 'T':
        #                 i = 1
        #     if self.sequence[-2] == 'G': i += 4
        #     self.is_conventional = i

        self.best_conv_var = self.is_conventional
        self.best_conv_obj = self
        if not self.variations:
            return
        for var in self.variations:
            var.conventional_version()
            if (var.is_conventional and conventional_class_rate(var.is_conventional)<conventional_class_rate(self.best_conv_var)) or (var.is_conventional and self.best_conv_var==0):
                self.best_conv_var = var.is_conventional
                self.best_conv_obj = var
        

    def nonconventional_version(self):
        def isR(N):
            if N in ["G","A"]: return True
            return False
        def isY(N):
            if N in ["C", "T"]: return True
            return False
        seq=self.sequence
        if self.next_exon and self.next_exon.sequence:
            nex=self.next_exon.sequence[:3]
            if len(nex)==2:
                nex=nex+" "
            elif len(nex)==1:
                nex=nex+"  "
        else: nex=None
        if self.prev_exon and self.prev_exon.sequence:
            prev = self.prev_exon.sequence[-1]
        else:
            prev = None
        i = 0
        try:
            if complimentary(seq[3],seq[-6]) and complimentary(seq[5],seq[-8]) and seq[4]=='A' and seq[-7]=='T': #9
                i=9
                if seq[3]=="C" and seq[-6]=='G': #8
                    i=8
                    if prev and nex and isY(prev) and isR(seq[0]) and isY(seq[-1]): #7/6
                        i = (isR(nex[0]) and 7) or (nex[2]=='C' and 6)
                        if isR(nex[0]) and nex[2]=='C': #3
                            i=3
                            if seq[5]=='G' and seq[-8]=='C': #2
                                i=2
                                if nex[1]=='A': #1
                                    i=1        
                    elif nex and isR(seq[0]) and isR(nex[0]) and nex[2]=='C': #5,4
                        if prev and isY(prev):
                            i=5
                        elif isY(seq[-1]):
                            i=4

            self.is_nonconventional = i
            self.best_nonconv_var = self.is_nonconventional
            self.best_nonconv_obj = self
            
            if not self.variations:
                return
            for var in self.variations:
                var.nonconventional_version()
                if (var.is_nonconventional and var.is_nonconventional < self.best_nonconv_var) or\
                        (var.is_nonconventional and self.best_nonconv_var == 0):
                    self.best_nonconv_var = var.is_nonconventional
        except Exception as exc:
            print(self.sequence)
            print(self.scaffold_name, self.scaffold_start, self.scaffold_end)
            print(self.prev_exon, self.next_exon)
            #raise exc

    def check_ML(self):
        if np.all(self.ML_characteristic == 0):
            self.calculate_ML_characteristic()
        # WAŻNE czy modele tak działają? czy z obu bierze się pozycje [0,1]?
        result_K = conventional_model.predict(self.ML_characteristic)[0]
        score_K = conventional_model.predict_proba(self.ML_characteristic)[0][1]
        result_NK = nonconventional_model.predict(self.ML_characteristic)[0]
        score_NK = nonconventional_model.predict_proba(self.ML_characteristic)[0][1]

        self.ML_class = [result_K, result_NK]
        self.ML_conv_score = score_K
        self.ML_nonconv_score = score_NK
        # if result_K == 1:
        #     self.ML_is_nonconventional = True
        # self.ML_nonconv_proba = proba
        # self.ML_best_conv_version = self
        # self.ML_best_nonconv_version = self
        for var in self.variations:
            var.check_ML()
        # for var in self.variations:
        #     c, p = var.check_ML()
        #     if p > self.ML_best_nonconv_version.ML_nonconv_proba:
        #         self.ML_best_nonconv_version = var
        #     if p < self.ML_best_conv_version.ML_nonconv_proba:
        #         self.ML_best_conv_version = var
        # if var_probas:
        #    if np.max(var_probas) > self.
        #    best_conv_v_ind, best_nonconv_v_ind = np.argmin(var_probas), np.argmax(var_probas)
        #    print(len(var_probas), (best_conv_v_ind, best_nonconv_v_ind))
        #    self.ML_best_conv_version = self.variations[best_conv_v_ind]
        #    self.ML_best_nonconv_version = self.variations[best_nonconv_v_ind]
        # return result, proba
    
    def is_ML_conv(self):
        """Returns whether the intron object is classified as conventional
        by the ML model"""
        return bool(self.ML_class[0])
    
    def is_ML_nonconv(self):
        """Returns whether the intron object is classified as conventional
        by the ML model"""
        return bool(self.ML_class[1])
    
    def get_ML_classification(self, k_or_nk):
        if k_or_nk == "K": return self.is_ML_conv()
        elif k_or_nk == "NK": return self.is_ML_nonconv()
        else:
            print("Only allowed types to check are K for conventional or NK for nonconventional")

    def calculate_ML_characteristic(self):
        if not self.prev_exon or not self.next_exon:
            return
        self.ML_characteristic = compute_intron_characteristics(self).reshape(1, -1)

    def add_manual_annotation(self, man_annotation, start, end):
        self.man_annotation = man_annotation
        if start == self.scaffold_start and end == self.scaffold_end:
            self.man_variant = self
        else:
            for var in self.variations:
                if start == var.scaffold_start and end == var.scaffold_end:
                    self.man_variant = var
                    break
        if not self.man_variant:
            print(self)
            print(man_annotation)

    def calculate_deka_score(self):
        penta_b = self.sequence[3:13]
        penta_e = self.sequence[-15:-5]
        return calculate_pairing(penta_b, penta_e)

    def add_test_annotation(self):
        # self.set_test_score()
        # self.test_score_max = self.test_score
        # self.test_score_min = self.test_score
        for var in self.variations:
            self.best_nonconv_obj = var
            var.set_test_score()
            if var.canonical_borders:
                self.test_k = True
            if var.test_score < self.test_score_min:
                self.test_score_min = var.test_score
                self.test_best_nk_var = var
        if self.test_best_nk_var.conserved_pairing_score > 5:
            self.test_nk = True
        if self.test_k:
            if self.test_nk:
                self.test_global_annotation = 'intron_I'
            else:
                self.test_global_annotation = 'intron_K'
        elif self.test_nk:
            self.test_global_annotation = 'intron_NK'
        else:
            self.test_global_annotation = 'intron_NN'
        # if self.test_score_max > 5 and self.test_score_min >= -5:
        #     self.test_global_annotation = 'intron_K'
        # elif self.test_score_min < -5 and self.test_score_max <= 5:
        #     self.test_global_annotation = 'intron_NK'
        # elif self.test_score_min < -5 and self.test_score_max > 5:
        #     if self.test_best_k_var.polypyrimidine_tract >= 0.6:
        #         if self.test_best_nk_var.conserved_pairing_score < 8:
        #             self.test_global_annotation = 'intron_K'
        #         else:
        #             self.test_global_annotation = 'intron_I'
        #     else:
        #         self.test_global_annotation = 'intron_I'
        # else:
        #     self.test_global_annotation = 'intron_NN'

    def set_test_score(self):
        if len(self.sequence) < 30:
            self.test_annotation = 'intron_NN'
        else:
            # self.test_score = 0
            # self.polypyrimidine_tract = calculate_pyrimidine_content(self.sequence[-12:-2])
            self.conserved_pairing_score = calculate_pairing(self.sequence[3:13], self.sequence[-15:-5])
            if self.sequence[:2] in ['GT', 'GC'] and self.sequence[-2:] == 'AG':
                self.canonical_borders = True
            else:
                self.canonical_borders = False
            # if self.canonical_borders:
            #     self.test_score += 10
            # if self.polypyrimidine_tract >= 0.6:
            #     self.test_score += 3
            # self.test_score -= self.conserved_pairing_score
            if self.canonical_borders:
                self.test_annotation = 'intron_K'
            elif self.conserved_pairing_score > 5:
                self.test_annotation = 'intron_NK'
            else:
                self.test_annotation = 'intron_NN'
    def is_man_K(self):
        return self.man_annotation == 'intron_K'
    
    def is_man_NK(self):
        return self.man_annotation == 'intron_NK'
    
    
    def noncanonical_structural_versions(self):
        s = self.sequence
        seq = self.prev_exon.sequence[-5:].lower() + s[:25]+'AAAAAAAAAA'+s[-25:] + self.next_exon.sequence[:5].lower()
        with open('seq_RNAfold.fasta', 'w') as f:
            f.write(">intron\n"+seq)
        #print(len(s), seq)
        stream = os.popen('RNAfold --noPS --auto-id --command=constr.txt < seq_RNAfold.fasta')
        output = stream.readlines()[2].split()[0]

        self.is_struct_nonconv_1, self.is_struct_nonconv_2 = False, False
        if not '(' in output[-25:] and not ')' in output[:25]:
            if output.count('(', 0, 25)>5:
                self.is_struct_nonconv_2 = True
                if self.next_exon.sequence[0] in ["G", "A"]:
                    self.is_struct_nonconv_1 = True

        for variation in self.variations:
            variation.noncanonical_structural_versions()

        if self.variations:
            self.variations_struct_1=sum([var.is_struct_nonconv_1 for var in self.variations])
            self.variations_struct_2=sum([var.is_struct_nonconv_2 for var in self.variations])


class Exon(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='', prev_exon=None, next_exon=None, prev_intron=None, next_intron=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.prev_exon = prev_exon
        self.next_exon = next_exon
        self.prev_intron = prev_intron
        self.next_intron = next_intron


class Transcript():
    def __init__(self, scaffold_name, start, end, sequence='', strand=''):
        self.scaffold_name = scaffold_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.strand = strand


def process_file(file_path):
    try:
        with open(file_path) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                yield [int(x) if x.isnumeric() else x for x in line.split()]
    except (IOError, OSError) as exc:
        print("Error opening / processing file")
        raise exc


def read_genome(file):
    genome = defaultdict(str)
    with open(file) as f:
        for line in f.readlines():
            if line[0] == '>':
                gene = line.strip()[1:]
            else:
                genome[gene] += line.strip()
    return genome
    

def read_genes(file):
    genes = {}  # slownik genow
    gene, exon = None, None
    prev = None
    for line in process_file(file):
        if line[0] == '#':
            continue
        if line[2] == 'transcript':
            if gene:
                genes[gene.name] = gene
            gene = Gene(line[0], line[3] - 1, line[4], name=line[11].strip('";'), strand=line[6], exons=[])
        elif line[2] == 'exon':
            exon = Exon(line[0], line[3] - 1, line[4], strand=line[6], prev_exon = prev)
            gene.append_exons(exon)
            if prev: prev.next_exon=exon
            prev=exon
    if gene:
        genes[gene.name] = gene
    return genes


def complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    letters = [complement_dict[base] for base in seq]
    return ''.join(letters)


def reverse_complement(seq):
    return complement(seq[::-1])

def complimentary(n1, n2):
    if {n1, n2} in [{'A', 'T'}, {'C', 'G'}, {'G', 'T'}]:
        return True
    else:
        return False

def conventional_class_rate(wersja):
    wagi={0:0, 1:1, 2:3, 3:5, 4:7, 5:2, 6:4, 7:6, 8:8}
    return wagi.get(wersja)

def create_genes(genome,genes):
    genome_eug = read_genome(genome)
    genes_eug = read_genes(genes)
    for name, gene in list(genes_eug.items()):
        gene.extract_sequence(genome_eug)
        gene.create_introns()
    return genome_eug,genes_eug

def liczenie_intronow(name, genome_eug,genes_eug, czy_fasta_gtf=False, czy_wykresy=False):
    if czy_fasta_gtf==True:
        genome,genes=create_genes(genome_eug, genes_eug)
        genome_eug=genome
        genes_eug=genes
    conv_count=0
    nconv_count=0
    both_count=0
    all_count=0
    non_count=0
    conventional_classes=[0,0,0,0,0,0,0,0]
    nonconventional_classes=[0,0,0,0,0,0,0,0,0]

    for genename, gene in list(genes_eug.items()):
        for intron in gene.introns:
            introns_to_assess.append(intron)
            their_characteristics.append(compute_intron_characteristics(intron))
            for var in intron.variations:
                introns_to_assess.append(var)
                their_characteristics.append(compute_intron_characteristics(var))
    print(len(introns_to_assess), len(their_characteristics))
    print("Finished predicting ML classes for introns!")
    assert len(their_characteristics) > 0
    
    predictions_conv = conventional_model.predict(their_characteristics)
    predictions_nonconv = nonconventional_model.predict(their_characteristics)
    probas_conv = conventional_model.predict_proba(their_characteristics)
    probas_nonconv = nonconventional_model.predict_proba(their_characteristics)
    #@TODO tu jest coś nie tak, gdzieś zwracane jest to samo dla obu model
    
    #probas = loaded_model.decision_function(their_characteristics)
    for (i, pred_nc, pred_c, prob_c, prob_nc) in\
        zip(introns_to_assess, predictions_nonconv, predictions_conv, probas_conv, probas_nonconv):
        i.ML_nonconv_score = prob_nc[1]
        i.ML_conv_score = prob_c[0]
        i.ML_class = (pred_c, pred_nc)
    #for intron in introns_to_assess:
        var_probas_conv = [v.ML_conv_score for v in i.variations]
        var_probas_nonconv = [v.ML_nonconv_score for v in i.variations]
        i.ML_best_conv_version = i.variations[np.argmin(var_probas_conv)] if len(var_probas_conv) else i
        i.ML_best_nonconv_version = i.variations[np.argmax(var_probas_nonconv)] if len(var_probas_nonconv) else i
    return

'''
if np.all(self.ML_characteristic == 0):
            self.calculate_ML_characteristic()
        result_K = conventional_model.predict(self.ML_characteristic)[0]
        score_K = conventional_model.predict_proba(self.ML_characteristic)[0][1]
        result_NK = nonconventional_model.predict(self.ML_characteristic)[0]
        score_NK = nonconventional_model.predict_proba(self.ML_characteristic)[0][1]

        self.ML_class = [result_K, result_NK]
        self.ML_conv_score = score_K
        self.ML_nonconv_score = score_NK
        # if result_K == 1:
        #     self.ML_is_nonconventional = True
        # self.ML_nonconv_proba = proba
        # self.ML_best_conv_version = self
        # self.ML_best_nonconv_version = self
        for var in self.variations:
            var.check_ML()
        # for var in self.variations:
        #     c, p = var.check_ML()
        #     if p > self.ML_best_nonconv_version.ML_nonconv_proba:
        #         self.ML_best_nonconv_version = var
        #     if p < self.ML_best_conv_version.ML_nonconv_proba:
        #         self.ML_best_conv_version = var
        # if var_probas:
        #    if np.max(var_probas) > self.
        #    best_conv_v_ind, best_nonconv_v_ind = np.argmin(var_probas), np.argmax(var_probas)
        #    print(len(var_probas), (best_conv_v_ind, best_nonconv_v_ind))
        #    self.ML_best_conv_version = self.variations[best_conv_v_ind]
        #    self.ML_best_nonconv_version = self.variations[best_nonconv_v_ind]
        # return result, proba
        # '''

def weigh_pairings(nn1, nn2):
    nn = nn1+nn2
    if nn in ["AT", "TA"]: return 0.5
    elif nn in ["GC", "CG"]: return 1
    elif nn in ["GT", "TG"]: return 0.375
    else: return 0

def get_introns_demulti(genes):
    """
    Get lists of conventional and nonconventional introns without repetitions
    from manual annotations.

    Args:
        genes (_type_): _description_

    Returns:
        introns_K_demulti (_type_): a list of unique conventional introns
        introns_NK_demulti (_type_): a list of unique nonconventional introns
    """
    introns_NK, introns_K, introns_I = [], [], []

    for _, gene_obj in genes.items():
        for intron in gene_obj.introns:
            if intron.man_annotation == 'intron_NK':
                introns_NK.append(intron)
            if intron.man_annotation == 'intron_K':
                introns_K.append(intron)
            if intron.man_annotation == 'intron_I':
                introns_I.append(intron)
    print('Intron_K', len(introns_K))
    print('Intron_NK', len(introns_NK))
    print('intron_I', len(introns_I))
    
    def demultiply_group(introns_X):
        introns_X_demulti = []
        single, double = [], {}
        for intron in introns_X:
            if intron.sequence in single:
                for i in range(len(single)):
                    if intron.sequence == single[i]:
                        name_1, name_2 = str(introns_X_demulti[i]), str(intron)
                        pos_1, pos_2 = name_1.split()[1:], name_2.split()[1:]
                        if pos_1 != pos_2:
                            if name_1 not in double.keys():
                                double[name_1] = [name_1, name_2]
                            else:
                                double[name_1].append(name_2)
            else:
                single.append(intron.sequence)
                introns_X_demulti.append(intron)
        return introns_X_demulti
    
    introns_K_demulti = demultiply_group(introns_K)
    introns_NK_demulti = demultiply_group(introns_NK)
                
    print('introns_K_demulti', len(introns_K_demulti))
    print('introns_NK_demulti', len(introns_NK_demulti))

    return introns_K_demulti, introns_NK_demulti

def conventional_class_rate(wersja):
    wagi={0:0, 1:1, 2:3, 3:5, 4:7, 5:2, 6:4, 7:6, 8:8}
    return wagi.get(wersja)

def create_genes(genome,genes):
    genome_eug = read_genome(genome)
    genes_eug = read_genes(genes)
    for name, gene in list(genes_eug.items()):
        gene.extract_sequence(genome_eug)
        gene.create_introns()
    return genome_eug,genes_eug

def liczenie_intronow(name, genome_eug,genes_eug, czy_fasta_gtf=False, czy_wykresy=False):
    if czy_fasta_gtf==True:
        genome,genes=create_genes(genome_eug, genes_eug)
        genome_eug=genome
        genes_eug=genes
    conv_count=0
    nconv_count=0
    both_count=0
    all_count=0
    non_count=0
    conventional_classes=[0,0,0,0,0,0,0,0]
    nonconventional_classes=[0,0,0,0,0,0,0,0,0]

    for genename, gene in list(genes_eug.items()):
        for intron in gene.introns:
            all_count+=1

            if intron.best_nonconv_var and not intron.best_conv_var:
                nonconventional_classes[intron.best_nonconv_var-1]+=1
                nconv_count+=1
            elif intron.best_conv_var and not intron.best_nonconv_var:
                conventional_classes[intron.best_conv_var-1]+=1
                conv_count+=1
            elif intron.best_conv_var and intron.best_nonconv_var:
                both_count+=1
            else:
                non_count+=1
    print("\nKonwencjonalne: %d, niekonwencjonalne: %d, oba: %d inne: %d, wszystkie: %d" %(conv_count,nconv_count, both_count,non_count,all_count))
    stats=[name, conventional_classes, nonconventional_classes, conv_count,nconv_count, both_count,non_count,all_count]
    #if czy_wykresy:
    #    wykresy_liczenia(stats)
    return stats

def wykresy_liczenia(stats, wykresy=("K"), tabele=True, save_wykresy=False):
    if wykresy:
        colors=[np.random.rand(3,) for i in range(len(stats))]
        for i in wykresy:
            if i not in ["K", "NK", "KNBN"]: print("Złe wartości wykresów")
    counts, klasy = [], []
    for stat in stats:
        name, conventional_classes, nonconventional_classes, conv_count,nconv_count, both_count,non_count,all_count=stat
        klasy.append([name,[i/all_count for i in conventional_classes],[i/all_count for i in nonconventional_classes]])
        counts.append([name, conv_count/all_count,nconv_count/all_count, both_count/all_count,non_count/all_count,all_count])
    #print(klasy, counts)
    if save_wykresy:
        print(counts)
        print(klasy)
        return
    
    if "K" in wykresy:
        lista1=[i[1] for i in klasy]
        data=[[j*100 for j in l] for l in lista1]
        legenda=[i[0] for i in klasy]
        X = np.arange(8)
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_axes([0,0,1,1])
        width = 1/(len(data)+1)
        for r in range(len(data)):
            ax.bar(X + width*r, data[r], color = colors[r], width = width)
        ax.set_title('Klasy intronów konwencjonalnych')
        ax.set_xlabel('Klasy konwencjonalne')
        ax.set_ylabel('% intronów')
        ax.set_xticks(X)
        ax.set_xticklabels(('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4', 'Klasa 5', 'Klasa 6', 'Klasa 7', 'Klasa 8'))
        ax.legend(legenda)
        if save_wykresy: plt.savefig("klasy_konwencjonalne.png")
    
    if "NK" in wykresy:
        lista1=[i[2] for i in klasy]
        data=[[j*100 for j in l] for l in lista1]
        legenda=[i[0] for i in klasy]
        Y = np.arange(9)
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_axes([0,0,1,1])
        width = 1/(len(data)+1)
        for r in range(len(data)):
            ax.bar(Y + width*r, data[r], color = colors[r], width = width)
        ax.set_title('Klasy intronów niekonwencjonalnych')
        ax.set_xlabel('Klasy niekonwencjonalne')
        ax.set_ylabel('% intronów')
        ax.set_xticks(Y)
        ax.set_xticklabels(('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4', 'Klasa 5', 'Klasa 6', 'Klasa 7', 'Klasa 8', 'Klasa 9'))
        ax.legend(legenda)
        if save_wykresy: plt.savefig("klasy_niekonwencjonalne.png")
    
    if "KNBN" in wykresy:
        data=[[j*100 for j in l[1:-1]] for l in counts]
        #print(data)
        labels=[i[0] for i in counts]
        Z = np.arange(4)
        fig = plt.figure(figsize=(7,5))
        ax = fig.add_axes([0,0,1,1])
        width = 1/(len(data)+1)
        for r in range(len(data)):
            ax.bar(Z + width*r, data[r], color = colors[r], width = width)
        ax.set_title('Ilości wszystkich intronów')
        ax.set_xlabel('Przynależność do klas')
        ax.set_ylabel('% intronów')
        ax.set_xticks(Z)
        ax.set_xticklabels(('Konwencjonalne','Niekonwencjonalne', 'Oba','Żadne'))
        ax.legend(labels)
    if tabele: print(tabulate(stats, headers=['nazwa', 'klasy konw.', 'klasy niekonw.', 'konwencjonalne', 'niekonwencjonalne', 'oba', 'żadne', "wszystkie"]))
    
    return
    
    