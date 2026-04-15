import os
from copy import copy
import numpy as np
from src.tools import auto_attr_check, complimentary, calculate_pairing, \
    calculate_gene_ends_from_scaffold, compute_intron_characteristics, \
        calculate_pyrimidine_content, conventional_class_rate
from src.models import nonconventional_model, conventional_model
from .genomic_sequence import GenomicSequence
from .gene import Gene

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
    :param is_conventional: (int) Optional, number of the conventional class the intron belongs to;
    if isn't conventional then 0.
    :param is_nonconventional: (int) Optional, number of the nonconventional class the intron belongs to;
    if isn't nonconventional then 0.
    :param best_conv_var: (int) Optional, unique to main intron, absent in var
    iations; number of the best conventional class out of all variants.
    :param best_nonconv_var: (int) Optional, unique to main intron, absent in variants;
    number of the best nonconventional class out of all variants.
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
    variants = list
    is_conventional = int
    is_nonconventional = int
    variations_struct_1 = int
    variations_struct_2 = int
    is_struct_nonconv_1 = int
    is_struct_nonconv_2 = int
    man_annotation = str
    test_annotation = str

    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None,
                 gene=None, support=None, margin_left=0, margin_right=0,
                 margin_left_seq='', margin_right_seq='', prev_exon=None, next_exon=None,
                 man_annotation='', test_annotation='', test_global_annotation='',
                 mother_of_intron=None, test_score=0, gc_content=0.0, polypyrimidine_tract=0.0):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end,
                                 sequence=sequence, strand=strand)
        self.gene = gene
        self.support = support
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_left_seq = margin_left_seq
        self.margin_right_seq = margin_right_seq
        self.variants = []
        self.mother_of_intron = mother_of_intron
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
        self.is_struct_nonconv_1 = False
        self.is_struct_nonconv_2 = False
        self.variants_struct_1 = None
        self.variants_struct_2 = None
        # TODO przy zmienianiu podstawowego intronu trzeba
        # przepisac best_(non)conv_var i liste wariacji
        #      - warianty ich nie maja i zawsze maja nie miec
        if self.strand == '-':
            self.gene_start, self.gene_end = -self.scaffold_end + self.gene.scaffold_end,\
                                             -self.scaffold_start + self.gene.scaffold_end
        elif self.gene and self.strand == '+':
            self.gene_start, self.gene_end = self.scaffold_start - self.gene.scaffold_start,\
                                             self.scaffold_end - self.gene.scaffold_start
        self.gene_start, self.gene_end = calculate_gene_ends_from_scaffold(self.gene.scaffold_start, self.gene.scaffold_end,
                               self.scaffold_start, self.scaffold_end, self.strand)
        if not self.scaffold_start < self.scaffold_end:
            raise ValueError(f"[INTRON] self.scaffold_start={self.scaffold_start}, self.scaffold_end={self.scaffold_end}")
        if (self.strand=='+' and not self.gene_start < self.gene_end) or \
            (self.strand=='-' and not self.gene_start > self.gene_end):
            raise ValueError(f"[INTRON] strand is {self.strand}, intron.gene_start={self.gene_start}, intron.gene_end={self.gene_end}")
        if not self.scaffold_name == self.gene.scaffold_name:
            raise ValueError(f"[INTRON] Wrong scaffonds, gene is {self.gene.scaffold_name}, intron is {self.scaffold_name}")

    def __str__(self):
        return f"""
    Intron in gene: {self.scaffold_name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}
    gene loc: {self.gene_start}-{self.gene_end}"""
    
    
    def movable_boundary_no_margins(self):
        # moglem zepsuc
        # ja zaraz jeszcze bardziej popsuje
        """
        Check if there are repeats on the intron junctions
        so the intron position could be shifted without changing
        transcript sequence. If there are, add new possible introns to self.variants.
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
                new_variation = Intron(self.scaffold_name, gene=self.gene, strand=self.strand, \
                    scaffold_start=self.scaffold_start - i, scaffold_end=self.scaffold_end - i, \
                    sequence=new_seq, prev_exon=new_prev_exon, next_exon=new_next_exon)

            else: #checking to the right
                if i + 1 > len(mrs) or i + 1 > len(self.sequence) or self.sequence[i] != mrs[i]:
                    break
                new_seq=self.sequence[i:]+mrs[:i]
                new_prev_exon.sequence = new_prev_exon.sequence+self.sequence[:i]
                new_prev_exon.scaffold_end = new_prev_exon.scaffold_end+i
                new_next_exon.sequence = new_next_exon.sequence[i:]
                new_next_exon.scaffold_start = new_next_exon.scaffold_start+i
                #creating new variation moved to the right
                new_variation = Intron(self.scaffold_name, gene=self.gene, strand=self.strand, \
                    scaffold_start=self.scaffold_start + i, scaffold_end=self.scaffold_end + i, \
                    sequence=new_seq, prev_exon=new_prev_exon, next_exon=new_next_exon)
            self.variants.append(new_variation)
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
        secondary structure in specific positions. Also check if variants with shifted junctions may be
        nonconventional."""

        left_anchor = self.sequence[3: 5]
        right_anchor = self.sequence[-7:-5]

        # try:
        if complimentary(left_anchor[0], right_anchor[1]) and complimentary(left_anchor[1], right_anchor[0]):
            return True
        # checking variants
        if len(self.variants) > 0:
            for son in self.variants:
                if son.check_nonconventional():
                    return True
        return False

    def conventional_version(self):
        if self.sequence[0:2] in ['GT', 'GC'] and self.sequence[-2:] == 'AG':
            i = 4
        else:
            return
        
        if not self.prev_exon or not self.next_exon:
            raise ValueError("why does the intron have no neighbouring exons?")
        mls, mrs = self.prev_exon.sequence[-3:], self.next_exon.sequence[:3]
        
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

        self.best_conv_var = self.is_conventional
        self.best_conv_obj = self
        if not self.variants:
            return
        for var in self.variants:
            var.conventional_version()
            if (var.is_conventional and conventional_class_rate(var.is_conventional)<conventional_class_rate(self.best_conv_var)) \
                or (var.is_conventional and self.best_conv_var==0):
                self.best_conv_var = var.is_conventional
                self.best_conv_obj = var
        
    # def nonconventional_version(self):
    #     def isR(N):
    #         if N in ["G","A"]:
    #             return True
    #         return False

    #     def isY(N):
    #         if N in ["C", "T"]:
    #             return True
    #         return False

    #     seq = self.sequence
    #     if self.next_exon and self.next_exon.sequence:
    #         nex = self.next_exon.sequence[:3]
    #         nex = nex + ' ' * (3 - len(nex))
    #     else:
    #         nex = None
    #     if self.prev_exon and self.prev_exon.sequence:
    #         prev = self.prev_exon.sequence[-1]
    #     else:
    #         prev = None
    #     i = 0
    #     try:
    #         if complimentary(seq[3], seq[-6]) and complimentary(seq[4], seq[-7]):
    #             i = 11
    #             if complimentary(seq[5], seq[-8]):  # 10
    #                 i = 10
    #                 if seq[4] == 'A' and seq[-7] == 'T':  # 9
    #                     i = 9
    #                     if seq[3] == "C" and seq[-6] == 'G':  # 8
    #                         i = 8
    #                         if prev and nex and isY(prev) and isR(seq[0]) and isY(seq[-1]):  # 7/6
    #                             i = (isR(nex[0]) and 7) or (nex[2] == 'C' and 6)
    #                             if isR(nex[0]) and nex[2] == 'C':  # 3
    #                                 i = 3
    #                                 if seq[5] == 'G' and seq[-8] == 'C':  # 2
    #                                     i = 2
    #                                     if nex[1] == 'A':  # 1
    #                                         i = 1
    #                         elif nex and isR(seq[0]) and isR(nex[0]) and nex[2] == 'C':  # 5,4
    #                             if prev and isY(prev):
    #                                 i = 5
    #                             elif isY(seq[-1]):
    #                                 i = 4

    #         self.is_nonconventional = i
    #         self.best_nonconv_var = self.is_nonconventional
    #         if not self.variants:
    #             return
    #         for var in self.variants:
    #             var.nonconventional_version()
    #             if (var.is_nonconventional and var.is_nonconventional < self.best_nonconv_var) or\
    #                     (var.is_nonconventional and self.best_nonconv_var == 0):
    #                 self.best_nonconv_var = var.is_nonconventional
    #     except Exception as exc:
    #         print(self.sequence)
    #         print(self.scaffold_name, self.scaffold_start, self.scaffold_end)
    #         print(self.prev_exon, self.next_exon)
    #         #raise exc
    
    def nonconventional_version(self):
        def isR(N):
            if N in ["G","A"]:
                return True
            return False

        def isY(N):
            if N in ["C", "T"]:
                return True
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
            if complimentary(seq[3],seq[-6]) and complimentary(seq[5],seq[-8]) \
                and seq[4]=='A' and seq[-7]=='T': #9
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

            if not self.variants:
                return
            for var in self.variants:
                var.nonconventional_version()
                if (var.is_nonconventional and var.is_nonconventional < self.best_nonconv_var) or\
                        (var.is_nonconventional and self.best_nonconv_var == 0):
                    self.best_nonconv_var = var.is_nonconventional
        except Exception as exc:
            print(self.sequence)
            print(self.scaffold_name, self.scaffold_start, self.scaffold_end)
            print(self.prev_exon, self.next_exon)
            print(exc)
            #pass

    def check_ML(self):
        if np.all(self.ML_characteristic == 0):
            self.calculate_ML_characteristic()
        # WAZNE czy modele tak dzialaja? czy z obu bierze sie pozycje [0,1]?
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
        for var in self.variants:
            var.check_ML()
        # for var in self.variants:
        #     c, p = var.check_ML()
        #     if p > self.ML_best_nonconv_version.ML_nonconv_proba:
        #         self.ML_best_nonconv_version = var
        #     if p < self.ML_best_conv_version.ML_nonconv_proba:
        #         self.ML_best_conv_version = var
        # if var_probas:
        #    if np.max(var_probas) > self.
        #    best_conv_v_ind, best_nonconv_v_ind = np.argmin(var_probas), np.argmax(var_probas)
        #    print(len(var_probas), (best_conv_v_ind, best_nonconv_v_ind))
        #    self.ML_best_conv_version = self.variants[best_conv_v_ind]
        #    self.ML_best_nonconv_version = self.variants[best_nonconv_v_ind]
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
        if k_or_nk == "K":
            return self.is_ML_conv()
        elif k_or_nk == "NK":
            return self.is_ML_nonconv()
        print("Only allowed types to check are K for conventional or NK for nonconventional")

    def calculate_ML_characteristic(self):
        if not self.prev_exon or not self.next_exon:
            return
        self.ML_characteristic = compute_intron_characteristics(self).reshape(1, -1)

    def add_manual_annotation(self, man_annotation, start, end):
        '''To be used when creating manual introns from gb file.'''
        self.man_annotation = man_annotation
        if start == self.scaffold_start and end == self.scaffold_end:
            self.man_variant = self
        else:
            for var in self.variants:
                if start == var.scaffold_start and end == var.scaffold_end:
                    self.man_variant = var
                    break
        if not self.man_variant:
            print("[add_manual_annotation] No manual variants", self)
            print(man_annotation)

    def calculate_deka_score(self):
        '''calculate_pairing btwn ends 3:13 and -15:-5'''
        penta_b = self.sequence[3:13]
        penta_e = self.sequence[-15:-5]
        return calculate_pairing(penta_b, penta_e)

    def add_test_annotation(self):
        '''
        Sets the following 7 attributes:
        - test_k - True if any variant has canonical borders
        - test_nk - True if the "least canonical" variant has high (>5) conserved_pairing_score
        - test_score_min/max - keeps track of the lowest/highest test_score btwn all variants
        - test_best_nk_var - variant with the lowest test_score (lowest canonicality score?)
        - test_best_k_var - variant with the highest test_score
        - test_global_annotation - intron type for the whole intron incl. variants
        '''
        self.set_test_score()
        self.test_score_max = self.test_score
        self.test_score_min = self.test_score
        for var in self.variants:
            var.set_test_score()
            if var.canonical_borders:
                self.test_k = True
            if var.test_score < self.test_score_min:
                self.test_score_min = var.test_score
                self.test_best_nk_var = var
            elif var.test_score > self.test_score_max:
                self.test_score_max = var.test_score
                self.test_best_k_var = var
        if self.test_best_nk_var.conserved_pairing_score > 5:
            self.test_nk = True

        if self.test_k and self.test_nk:
            self.test_global_annotation = 'intron_I'
        elif self.test_k and not self.test_nk:
            self.test_global_annotation = 'intron_K'
        elif not self.test_k and self.test_nk:
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
        '''
        Sets the following attributes:
        - conserved_pairing_score - calculate_pairing btwn ends 3:13 and -15:-5
        - canonical_borders - G{T/C}...AG
        - test_annotation [intron_K, intron_NK, intron_NN]
            intron_K - when intron has canonical borders
            intron_NK - when intron has enough pairing btwn ends (conserved_pairing_score>5)
            intron_NN - if intron is shorter than 30 and in all other cases
        - test_score [int] - some score of how canonical the intron is??
        - polypyrimidine_tract [0-1] - how much of the end of the sequence (-12:-2) are pyrimidines
        '''
        if len(self.sequence) < 30:
            self.test_annotation = 'intron_NN'
            return

        self.test_score = 0
        self.polypyrimidine_tract = calculate_pyrimidine_content(self.sequence[-12:-2])
        self.conserved_pairing_score = self.calculate_deka_score()
        self.canonical_borders = (self.sequence[:2] in ['GT', 'GC'] and self.sequence[-2:]=='AG')
        if self.canonical_borders:
            self.test_score += 10
        if self.polypyrimidine_tract >= 0.6:
            self.test_score += 3
        self.test_score -= self.conserved_pairing_score

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
        with open('seq_RNAfold.fasta', 'w', encoding='utf-8') as f:
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

        for variation in self.variants:
            variation.noncanonical_structural_versions()

        if self.variants:
            self.variants_struct_1=sum(var.is_struct_nonconv_1 for var in self.variants)
            self.variants_struct_2=sum(var.is_struct_nonconv_2 for var in self.variants)
