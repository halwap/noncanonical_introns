import os
import pickle
from bisect import insort
import warnings
from re import search
from copy import copy
from collections import defaultdict
import os
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO
from tqdm import tqdm
from tabulate import tabulate
import Bio
from typing import TextIO
from itertools import pairwise

warnings.filterwarnings("ignore",
                        #message="divide by zero encountered in divide",
                        category=UserWarning)



def HELP_load_default_genome_genes(species=None, with_reversed=False, if_classifier_mixed=False):
    '''@TODO
    help, not to be used eventually'''

    path_a = "fastas_and_gff/Alignments/"
    if not species:
        genome, genes = create(path_a+'test_edited_alignments_all_clean_3.fasta',
                            path_a+'test_edited_alignments_all_clean_3.gff', 'manual')
        add_manual_annotation_from_gb(path_a+'alignments_all_clean_3_genes_only.gb', genes)
        predict_all_introns(genes, if_mixed=if_classifier_mixed)
        return genome, genes

    #path="/mnt/archive/euglena_genomy/"
    path="./fastas_and_gff/Alignments/"
    #path_rev="./fastas_and_gff/reversed_fasta_gtf/"
    path_rev="./fastas_and_gff/Alignments/"

    if species == "EG":
        genome_file=path+"gracilis_dbg2olc_nogaps.fasta"
        genes_file=path+"gracilis_stringtie_strand_informed_nogaps.gtf"
        genes_rev_file=path_rev+"gracilis_stringtie_strand_informed_reversed_nogaps.gtf"

    elif species == "EH":
        genome_file=path+'hiemalis_rascaf_nogaps.fasta'
        genes_file=path+'hiemalis_stringtie_strand_informed3_nogaps.gtf'
        genes_rev_file=path_rev+"hiemalis_stringtie_strand_informed3_reversed_nogaps.gtf"

    elif species == "EL":
        genome_file=path+'longa_rascaf_nogaps.fasta'
        genes_file=path+'longa_stringtie_strand_informed_nogaps.gtf'    
        genes_rev_file=path_rev+"longa_stringtie_strand_informed_reversed_nogaps.gtf"

    elif species == "bugtest":
        genome_file = './fastas_and_gff/bugtest.fasta'
        genes_file = './fastas_and_gff/bugtest.gtf'
        genes_rev_file = './fastas_and_gff/bugtest_reversed.gtf'

    else:
        print(f"Incorrect species input: {species}")

    if not with_reversed:
        r = create(genome_file, genes_file, "stringtie")
    else:
        r = create_both(genome_file, genes_file, genes_rev_file, "stringtie")

    return r #genome_eug, genes_eug, [genes_eug_rev]



def getter_setter_gen(name, type_):
    def getter(self):
        return getattr(self, "__" + name)

    def setter(self, value):
        if not isinstance(value, type_) and value is not None:
            raise TypeError(f"{name} attribute must be set to an instance of {type_}")
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
        if not self.scaffold_start<self.scaffold_end:
            raise ValueError(f"Incorrect start & end coordinates: start={self.scaffold_start}, end={self.scaffold_end} in scaffold {scaffold_name}")
        if sequence and len(sequence) != self.length():
            raise ValueError(f'Incorrect sequence length in {self}: len(seq) = {len(sequence)}, self.length={self.length()}')
        self.sequence = sequence
        if not strand or (strand and strand not in ['+', '-', '.']):
            raise ValueError(f'Strand can only be +, - or ., is {strand}')
        self.strand = strand

    def length(self):
        return self.scaffold_end - self.scaffold_start
    
    def __repr__(self):
        return ' '.join([self.scaffold_name, str(self.scaffold_start), str(self.scaffold_end)])
        
    def __str__(self):
        return self.__repr__()


class Gene(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='',
                 transcript=None, exons=None, introns=None,
                 working_exons=None, working_introns=None,
                 name='', coverage=0):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=sequence, strand=strand)
        self.transcript = transcript
        self.working_exons = working_exons if working_exons else []
        self.working_introns = working_introns if working_introns else []
        self.exons = [] if not exons else exons
        self.introns = [] if not introns else introns
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


    # def add_exons(self, genes_data_type):
    #     if self.strand == '+':
    #         self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
    #     elif self.strand == '-':
    #         self.working_exons.sort(key=lambda _exon: _exon.scaffold_start, reverse=True)
    #     else:
    #         raise ValueError('gene strand not in [+, -]'.format(self.strand))
    #     prev = None
    #     for exon in self.working_exons:
    #         if prev:
    #             if genes_data_type == 'manual':
    #                 if prev.scaffold_end == exon.scaffold_start:
    #                     prev.merge_exons(exon)
    #                     continue
    #             prev.next_exon = exon
    #             exon.prev_exon = prev
    #         self.exons.append(exon)
    #         prev = exon
    
    def __str__(self):
        return f"""
    Gene {self.name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}"""
    
    def append_exons(self, exon):
        '''Adds exon to the exons list in the correct order'''
        insort(self.exons, exon, key=lambda _exon: _exon.scaffold_start)
        self.exons.append(exon)


    def add_exons(self, genes_data_type):
        '''Adding exons to the gene.'''
        if len(self.working_exons)==0:
            return
        if len(self.working_exons)==1:
            self.exons = self.working_exons
            return
                
        #self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
        self.working_exons = sorted(self.working_exons, key=lambda _exon: _exon.scaffold_start)

        exon_gen = (e for e in self.working_exons)
        prev_exon = None
        try:
            curr_exon = next(exon_gen)
        except StopIteration:
            return

        def exon_s_e(exon):
            if exon is None:
                return '_:_'
            return f"{exon.scaffold_start}:{exon.scaffold_end}"
        
        while curr_exon:
            #print(f"[loop start] prev_exon {exon_s_e(prev_exon)} curr_exon {exon_s_e(curr_exon)}")
            if prev_exon:
                if genes_data_type == 'manual' and \
                    prev_exon.scaffold_end == curr_exon.scaffold_start:
                    #merging exons if there is no intron between them
                    new_exon = copy(prev_exon)
                    new_exon.sequence = prev_exon.sequence + curr_exon.sequence
                    new_exon.scaffold_end = curr_exon.scaffold_end
                    new_exon.next_exon = curr_exon.next_exon
                    #print(f"merging {exon_s_e(prev_exon)} with {exon_s_e(curr_exon)} into {exon_s_e(new_exon)}")
                    if len(self.exons)>0:
                        self.exons[-1] = new_exon
                    else:
                        self.exons.append(new_exon)
                    prev_exon, curr_exon = new_exon, new_exon.next_exon
                else:
                    prev_exon.next_exon = curr_exon
                    curr_exon.prev_exon = prev_exon
                    self.exons.append(curr_exon)
                    #print(f"appended {exon_s_e(curr_exon)}")
                    prev_exon = curr_exon
            else:
                self.exons.append(curr_exon)
                prev_exon = curr_exon
            try:
                curr_exon = next(exon_gen)
            except StopIteration:
                break
        
        if len(self.working_exons)>0 and len(self.exons) == 0:
            print("[add_exons] why were none of the following exons created??")
            for e in self.working_exons:
                print('\t', e)

        if genes_data_type != 'manual' and len(self.exons) != len(self.working_exons):
            raise ValueError(f"[add_exons] Created {len(self.exons)} exons, while should be {len(self.working_exons)}")
        
    
    def append_introns(self, intron):
        self.introns.append(intron)
    
    def create_intron_structures(self, path_working="./ViennaRNA-2.4.18/working_fastas/", path_RNAfold = '/usr/local/bin/'):
        
        p=plik.split('/')[-1].split('.')[0]
        command = [path_RNAfold+'RNAfold', '--noPS',  '--command=constr.txt', '-j8', '-i', path_working+p+'.fasta']
        with open(path_working+p+'.dbn', 'w') as output_file:
            _p = subprocess.run(command, stdout=output_file)


    # def add_exons(self, genes_data_type):
    #     if self.strand == '+':
    #         self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
    #     elif self.strand == '-':
    #         self.working_exons.sort(key=lambda _exon: _exon.scaffold_start, reverse=True)
    #     else:
    #         raise ValueError('gene strand not in [+, -]'.format(self.strand))
    #     prev = None
    #     for exon in self.working_exons:
    #         if prev:
    #             if genes_data_type == 'manual':
    #                 if prev.scaffold_end == exon.scaffold_start:
    #                     prev.merge_exons(exon)
    #                     continue
    #             prev.next_exon = exon
    #             exon.prev_exon = prev
    #         self.exons.append(exon)
    #         prev = exon
    
    def __str__(self):
        return f"""
    Gene {self.name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}"""
    
    def append_exons(self, exon):
        '''Adds exon to the exons list in the correct order'''
        insort(self.exons, exon, key=lambda _exon: _exon.scaffold_start)
        self.exons.append(exon)


    def add_exons(self, genes_data_type):
        '''Adding exons to the gene.'''
        if len(self.working_exons)==0:
            return
        if len(self.working_exons)==1:
            self.exons = self.working_exons
            return
                
        #self.working_exons.sort(key=lambda _exon: _exon.scaffold_start)
        self.working_exons = sorted(self.working_exons, key=lambda _exon: _exon.scaffold_start)

        exon_gen = (e for e in self.working_exons)
        prev_exon = None
        try:
            curr_exon = next(exon_gen)
        except StopIteration:
            return

        def exon_s_e(exon):
            if exon is None:
                return '_:_'
            return f"{exon.scaffold_start}:{exon.scaffold_end}"
        
        while curr_exon:
            #print(f"[loop start] prev_exon {exon_s_e(prev_exon)} curr_exon {exon_s_e(curr_exon)}")
            if prev_exon:
                if genes_data_type == 'manual' and \
                    prev_exon.scaffold_end == curr_exon.scaffold_start:
                    #merging exons if there is no intron between them
                    new_exon = copy(prev_exon)
                    new_exon.sequence = prev_exon.sequence + curr_exon.sequence
                    new_exon.scaffold_end = curr_exon.scaffold_end
                    new_exon.next_exon = curr_exon.next_exon
                    #print(f"merging {exon_s_e(prev_exon)} with {exon_s_e(curr_exon)} into {exon_s_e(new_exon)}")
                    if len(self.exons)>0:
                        self.exons[-1] = new_exon
                    else:
                        self.exons.append(new_exon)
                    prev_exon, curr_exon = new_exon, new_exon.next_exon
                else:
                    prev_exon.next_exon = curr_exon
                    curr_exon.prev_exon = prev_exon
                    self.exons.append(curr_exon)
                    #print(f"appended {exon_s_e(curr_exon)}")
                    prev_exon = curr_exon
            else:
                self.exons.append(curr_exon)
                prev_exon = curr_exon
            try:
                curr_exon = next(exon_gen)
            except StopIteration:
                break
        
        if len(self.working_exons)>0 and len(self.exons) == 0:
            print("[add_exons] why were none of the following exons created??")
            for e in self.working_exons:
                print('\t', e)

        if genes_data_type != 'manual' and len(self.exons) != len(self.working_exons):
            raise ValueError(f"[add_exons] Created {len(self.exons)} exons, while should be {len(self.working_exons)}")
        

    def extract_sequence(self, genome):
        # elif self.strand == '-':
        #     sequence = genome[self.scaffold_name][self.scaffold_end:self.scaffold_start]
        # else:
        #     raise Exception('cojest')
        scaffold_seq = genome[self.scaffold_name]
        sequence = scaffold_seq[self.scaffold_start:self.scaffold_end]
        expanded_sequence, expansion_left, expansion_right = self.get_expanded_sequence(scaffold_seq)
        if self.strand == '-':
            self.sequence = Bio.Seq.reverse_complement(sequence)
            self.expanded_sequence = Bio.Seq.reverse_complement(expanded_sequence)
            self.expansion_left = expansion_right
            self.expansion_right = expansion_left
        else:
            self.sequence = sequence
            self.expanded_sequence = expanded_sequence
            self.expansion_right = expansion_right
            self.expansion_left = expansion_left
        # for exon in self.exons:
        #     if self.strand == '+':
        #         exon.sequence = self.sequence[exon.scaffold_start - self.scaffold_start:
        #                                       exon.scaffold_end - self.scaffold_start]
        #     elif self.strand == '-':
        #         exon.sequence = self.sequence[- exon.scaffold_end + self.scaffold_end:
        #                                       - exon.scaffold_start + self.scaffold_end]
        #     else:
        #         raise ValueError('gene strand: {} not in [+, -]'.format(self.strand))
        # for exon in self.exons:
        #     if self.strand == '+':
        #         exon.sequence = self.sequence[exon.scaffold_start - self.scaffold_start:
        #                                       exon.scaffold_end - self.scaffold_start]
        #     elif self.strand == '-':
        #         exon.sequence = self.sequence[-exon.scaffold_end + self.scaffold_end:
        #                                       -exon.scaffold_start + self.scaffold_end]
        #     else:
        #         print("Sequence extraction failed\t", self.name, exon.scaffold_name, exon.scaffold_start, exon.scaffold_end)
        #         # raise Exception('co jest')
        for exon in self.exons:
            exon_gene_start, exon_gene_end = calculate_gene_ends_from_scaffold(self.scaffold_start, self.scaffold_end,
                                                                    exon.scaffold_start, exon.scaffold_end, self.strand)
            if self.strand=="-":
                exon.sequence = self.sequence[exon_gene_start:exon_gene_end:-1]
            else:
                exon.sequence = self.sequence[exon_gene_start:exon_gene_end]
                
        transcript_sequence = self.get_transcript_sequence()
        self.transcript = Transcript(self.scaffold_name, self.scaffold_start, self.scaffold_end,
                                     strand=self.strand, sequence=transcript_sequence)
        
        if not len(sequence) == abs(self.scaffold_start-self.scaffold_end):
            raise ValueError(f"""The extracted gene sequence is of wrong length:
                             len(seq)={len(sequence)}
                             start-end={self.scaffold_start-self.scaffold_end}""")

    def get_expanded_sequence(self, scaffold_seq):
        start = max(0, self.scaffold_start - 500)
        end = min(len(scaffold_seq), self.scaffold_end + 500)
        expanded_sequence = scaffold_seq[start:end]
        expansion_left = self.scaffold_start - start
        expansion_right = end - self.scaffold_end
        return expanded_sequence, expansion_left, expansion_right

    def get_transcript_sequence(self):
        # exons_seqs = []
        # for exon in self.exons:
        #     exons_seqs.append(exon.sequence)
        # sequence = ''.join([exon_seq for exon_seq in exons_seqs])
        # return sequence
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

    def predict_introns(self, if_mixed=False):
        '''Predicting introns for one gene.
        Only use with singular gene, if with whole genome - use predict_all_introns'''
        introns_to_assess, their_characteristics = [], []
        for intron in self.introns:
            introns_to_assess.append(intron)
            their_characteristics.append(compute_intron_characteristics(intron))
            for var in intron.variants:
                introns_to_assess.append(var)
                their_characteristics.append(compute_intron_characteristics(var))
        if len(their_characteristics) == 0:
            return
        
        if if_mixed:
            nonconv_predictions = get_nonconventional_model().predict(their_characteristics)
            for i, nc_pred in zip(introns_to_assess, nonconv_predictions):
                is_k = int(bool(i.is_conventional))
                i.ML_class = [is_k, nc_pred]
        else:
            # WAZNE nie pomylic indeksow, conv to 1 a nonconv to 0
            conv_predictions = 1-get_conventional_model().predict(their_characteristics)
            conv_scores = get_conventional_model().predict_proba(their_characteristics)[:, 1]
            nonconv_scores = get_nonconventional_model().predict_proba(their_characteristics)[:, 0]
            for i, c_pred, c_score, nc_pred, nc_score in \
                    zip(introns_to_assess, conv_predictions, conv_scores,
                        nonconv_predictions, nonconv_scores):
                i.ML_class = [c_pred, nc_pred]
                i.ML_conv_score = c_score
                i.ML_nonconv_score = nc_score
            for intron in self.introns:
                var_probas = [v.ML_nonconv_score for v in intron.variants]
                intron.ML_best_conv_version = intron.variants[np.argmin(var_probas)]\
                    if len(var_probas) else intron
                intron.ML_best_nonconv_version = intron.variants[np.argmax(var_probas)]\
                    if len(var_probas) else intron
        return

    # def create_introns(self):
    #     if len(self.exons) < 1:
    #         raise ValueError('No exons specified.')
    #     elif len(self.exons) > 1:
    #         self.introns = []
    #         start, end = 0, 0
    #         for exon in self.exons:
    #             prev_exon = exon.prev_exon
    #             if self.strand == '+':
    #                 end = exon.scaffold_start
    #             elif self.strand == '-':
    #                 start = exon.scaffold_end
    #             if start and end:
    #                 if self.strand == '+':
    #                     sequence = self.sequence[start - self.scaffold_start:end - self.scaffold_start]
    #                 elif self.strand == '-':
    #                     sequence = self.sequence[- end + self.scaffold_end: - start + self.scaffold_end]
    #                 intron = Intron(self.scaffold_name, scaffold_start=start, scaffold_end=end, strand=self.strand,
    #                                 sequence=sequence, gene=self, prev_exon=prev_exon, next_exon=exon)
    #                 self.introns.append(intron)
    #                 prev_exon.next_intron = intron
    #                 exon.prev_intron = intron
    #             if self.strand == '+':
    #                 start = exon.scaffold_end
    #             elif self.strand == '-':
    #                 end = exon.scaffold_start
    #         for intron in self.introns:
    #             intron.movable_boundary_no_margins()
    #             intron.conventional_version()
    #             intron.nonconventional_version()
    #             intron.add_test_annotation()

    def create_introns(self):
        '''Creates all introns within the gene.
        An intron's start is the end of an exon,
        and end is the start of a following exon.
        Before the first exon there is no intron.'''

        if len(self.exons) < 1:
            raise ValueError(f'[create_introns] No exons specified. {self}')
        if len(self.exons) == 1:
            #print("Only one exon in gene, no introns will be created")
            return

        self.introns = []
        intron_scaffold_start = self.exons[0].scaffold_end
        intron_scaffold_end = 0
        current_exon = self.exons[1]
        while current_exon:
            #creating intron that is btwn prev_exon and current_exon
            prev_exon = current_exon.prev_exon
            intron_scaffold_end = current_exon.scaffold_start
            intron_gene_start, intron_gene_end = calculate_gene_ends_from_scaffold(self.scaffold_start, \
                self.scaffold_end, intron_scaffold_start, intron_scaffold_end, self.strand)
            
            if self.strand=="-":
                intron_sequence = self.sequence[intron_gene_start:intron_gene_end:-1]
            else:
                intron_sequence = self.sequence[intron_gene_start:intron_gene_end]

            
            if len(intron_sequence)!=abs(intron_gene_start-intron_gene_end):
                raise ValueError(f"""gene sequence length: {len(self.sequence)}
                                 intron_gene_start:intron_gene_end {intron_gene_start}:{intron_gene_end}""")
            
            intr = Intron(self.scaffold_name, scaffold_start = intron_scaffold_start, \
                scaffold_end = intron_scaffold_end, \
                strand=self.strand, sequence=intron_sequence, gene=self,\
                prev_exon=prev_exon, next_exon=current_exon)
            #self.append_introns(intr)
            self.introns.append(intr)

            if prev_exon:
                prev_exon.next_intron = intr
            current_exon.prev_intron=intr

            prev_exon=current_exon
            # starting a new intron that follows current_exon,
            # which will be prev_exon in the following loop
            intron_scaffold_start = current_exon.scaffold_end

            current_exon = current_exon.next_exon
        
        for intron in self.introns:
            if not len(intron.sequence) == abs(intron.scaffold_start-intron.scaffold_end):
                raise ValueError(f"""The extracted intron sequence is of wrong length:
                    len(seq)={len(intron.sequence)}
                    start:end={abs(intron.scaffold_start-intron.scaffold_end)}, {abs(intron.gene_end-intron.gene_start)}""")
                
            if not (intron.prev_exon and intron.next_exon):
                print("[create_introns] created intron has no neighbouring exons", intron)
                print(self.exons)
            intron.movable_boundary_no_margins()
            intron.conventional_version()
            intron.nonconventional_version()
    
    def finalize_serialize(self, fd: TextIO):
        """
        For each intron within the gene, find the highest scored variant,
        and takes that variant to be the de facto intron.
        With the introns rectified, the entire gene is serialized to GFF format.
        The gene should have already scored introns.
        
        Parameters:
        fd -- The file to serialize to. A file descriptor, not a filename!
        """
        
        #This method is only for multi-exon genes; if the gene is single-exon,
        #use single_exon_serialize() instead
        if len(self.exons) == 1:
            self.single_exon_serialize(fd)
            return
        
        
        #Check whether the gene has introns at all
        if len(self.introns) == 0:
            raise ValueError(f"{self} does not contain introns")
        
        #Check whether the gene has the appropriate number of introns
        if len(self.introns) != len(self.exons)-1:
            raise ValueError(f"{self} contains incorrect number of introns")
        
        #Check whether introns have been scored
        #This check might be either overly or insufficiently exhaustive...
        if all( intron.ML_nonconv_score == 0 for intron in self.introns ):
            raise ValueError(f"Introns of {self} have not been scored")
        
        #Check whether the terminal exons' start & end match the gene's start
        #and end
        if self.exons[0].scaffold_start != self.scaffold_start or \
           self.exons[-1].scaffold_end != self.scaffold_end:
            raise ValueError(f"Exons of {self} unaligned with the gene")
        
        #Likewise for transcript
        if self.transcript.start != self.scaffold_start or \
           self.transcript.end != self.scaffold_end:
            raise ValueError(f"Transcript of {self} unaligned with the gene")
        
        
        
        #For each intron, find the highest scored variant
        #Will hold each intron's best variant
        best_var: list[Intron] = []
        
        
        for intron in self.introns:
            #Will maximalize score
            best_c_score: float = intron.ML_conv_score
            best_nc_score: float = intron.ML_nonconv_score
            #Index of best variant
            #-1 is the main intron, >=0 is an entry in `intron.variants'
            best_c_index: int = -1
            best_nc_index: int = -1
            
            for idx, variant in enumerate(intron.variants):
                if variant.ML_conv_score > best_c_score:
                    best_c_score = variant.ML_conv_score
                    best_c_index = idx
                if variant.ML_nonconv_score > best_nc_score:
                    best_nc_score = variant.ML_nonconv_score
                    best_nc_index = idx
            
            #Determine whether the absolute best score was a conventional one or
            #a nonconventional one
            #The conventional score takes priority over the nonconventional one
            #in case of ties (exceedingly unlikely, since the scores are floats)
            if best_c_score >= best_nc_score:
                #Index of variant with absolute highest variant
                best_index = best_c_index
            else:
                best_index = best_nc_index
            
            
            #Remember the best variant, and whether the conventional or
            #nonconventional score won
            best_var.append(intron if best_index == -1 else intron.variants[best_index])
        
        
        #With each intron's best variant found, it's time to serialize the gene
        #The important positions to keep track are the start/end of each intron,
        #plus the start and end of the first and last exon, respectively
        #Start of gene, and first exon
        start = self.scaffold_start + 1
        #End of gene, and last exon
        end = self.scaffold_end
        
        
        #Write gene feature
        self.emit_gff(fd, "gene", start, end)
        self.emit_gff(fd, "mRNA", start, end)
        
        
        #Get start and end end position of each intron
        intron_pos = [ (i.scaffold_start+1, i.scaffold_end) for i in best_var ]
        
        #Write initial exon
        self.emit_gff(fd, "exon", start, intron_pos[0][0]-1, 1)
        
        #Write medial exons
        for n, ((prev_start, prev_end), (next_start, next_end)) \
        in enumerate(pairwise(intron_pos), 2):
            self.emit_gff(fd, "exon", prev_end+1, next_start-1, n)
        
        #Write terminal exon
        n = len(self.exons)
        self.emit_gff(fd, "exon", intron_pos[-1][1]+1, end, len(self.exons))
        
        
        #Write introns
        for n in range(len(self.introns)):
            self.emit_gff(fd, "intron", intron_pos[n][0], intron_pos[n][1], n+1)
    
    
    def single_exon_serialize(self, fd: TextIO):
        """
        Serializes a single-exon gene to GFF format.
        
        Parameters:
        fd  The file to serialize to. A file descriptor, not a filename!
        """
        
        #Check whether the exon's start & end match the gene's start and end
        if self.exons[0].scaffold_start != self.scaffold_start or \
           self.exons[0].scaffold_end != self.scaffold_end:
            raise ValueError(f"Exon of {self} unaligned with the gene")
        
        #Likewise for transcript
        if self.transcript.start != self.scaffold_start or \
           self.transcript.end != self.scaffold_end:
            raise ValueError(f"Transcript of {self} unaligned with the gene")
        

        start = self.scaffold_start + 1
        #End of gene, and last exon
        end = self.scaffold_end
        

        #Write features
        self.emit_gff(fd, "gene", start, end)
        self.emit_gff(fd, "mRNA", start, end)
        self.emit_gff(fd, "exon", start, end)
    
    
    def emit_gff(self, fd: TextIO, type: str, start: int, end: int, nth: int=1):
        """
        Helper method for finalize_serialize() and single_exon_serialize();
        emits a single line of GFF describing a gene/transcript/exon/intron.
        
        Parameters:
        fd      TextIO  File to write to. A file descriptor, not a filename!
        type    str     Type of feature: "gene", "mRNA", "exon" or "intron"
        start   int     Start position
        end     int     End position
        nth     int     Ordinal number of exon/intron within gene
        """
        
        #Check feature type
        if type not in { "gene", "mRNA", "exon", "intron" }:
            raise ValueError(f"Invalid feature type to serialize: {type}")
        
        
        #Build attributes field
        attr: str = f"Name={self.name};ID={self.name}"
        
        #If this is a transcript, mark ID as such and set gene as parent
        if type == "mRNA":
            attr += f".mrna;Parent={self.name}"
        #If this is an exon/intron, mark ID as such and set transcript as parent
        elif type == "exon" or type == "intron":
            attr += f".mrna.{type}{nth};Parent={self.name}.mrna"
        
        
        #Build GFF record
        record: str = '\t'.join([
                self.scaffold_name,     #Scaffold
                ".",                    #Source of feature
                type,                   #Feature type
                str(start),             #Start position
                str(end),               #End position
                ".",                    #Score
                self.strand,             #Strand
                ".",                    #Phase
                attr                    #Attributes
            ]) + '\n'
        
        #Write to file
        fd.write(record)





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
        #Precalculate position bounds for the check
        left_bound = min(len(self.sequence), len(mls))
        right_bound = min(len(self.sequence), len(mrs)) - 1
        
        
        i = 1
        # start checking for repeats left from the junction
        check_left = True
        
        
        while True:
            
            if check_left: #checking to the left
                if i > left_bound or mls[-i] != self.sequence[-i]:
                    check_left = False
                    i = 0
                    continue

                new_prev_exon = copy(self.prev_exon)
                new_next_exon = copy(self.next_exon)
                
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
                if i > right_bound or self.sequence[i] != mrs[i]:
                    break
                
                new_prev_exon = copy(self.prev_exon)
                new_next_exon = copy(self.next_exon)
                
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
            return N=='G' or N=='A'
        
        def isY(N):
            return N=='C' or N=='T'
        
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
            #print(self.sequence)
            #print(self.scaffold_name, self.scaffold_start, self.scaffold_end)
            #print(self.prev_exon, self.next_exon)
            pass

    def check_ML(self):
        if np.all(self.ML_characteristic == 0):
            self.calculate_ML_characteristic()
        # WAZNE czy modele tak dzialaja? czy z obu bierze sie pozycje [0,1]?
        result_K = get_conventional_model().predict(self.ML_characteristic)[0]
        score_K = get_conventional_model().predict_proba(self.ML_characteristic)[0][1]
        result_NK = get_nonconventional_model().predict(self.ML_characteristic)[0]
        score_NK = get_nonconventional_model().predict_proba(self.ML_characteristic)[0][1]

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
        with open('seq_RNAfold.fasta', 'w') as f:
            fd.write(">intron\n"+seq)
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


class Exon(GenomicSequence):
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='',
                 gene=None, prev_exon=None, next_exon=None, prev_intron=None, next_intron=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end,
                                 sequence=sequence, strand=strand)
        if (prev_intron and self.scaffold_name!=prev_intron.scaffold_name) or \
            (next_intron and self.scaffold_name!=next_intron.scaffold_name):
            raise ValueError(f"""[EXON] Neighbouring introns are in different scaffolds:
                             prev i: {prev_intron}, current: {self}, next i: {next_intron}""")
        self.prev_intron = prev_intron
        self.next_intron = next_intron
        self.gene = gene
        if gene:
            self.check_gene()

        if (prev_exon and prev_exon.scaffold_name!=self.scaffold_name) \
            or (next_exon and next_exon.scaffold_name!=self.scaffold_name):
            raise ValueError(f"""[EXON] Neighbouring exons are in different scaffolds:
                             prev e: {prev_exon}, current: {self}, next e: {next_exon}""")
        if prev_exon and (prev_exon.scaffold_end!=scaffold_start+1):
            raise ValueError(f"[EXON] self.prev_exon.scaffold_end={prev_exon.scaffold_end} != self.scaffold_start+1={scaffold_start}+1")
        self.prev_exon = prev_exon
        self.next_exon = next_exon

        if self.next_exon and (self.scaffold_end+1!=self.next_exon.scaffold_start):
            raise ValueError(f"[EXON] self.scaffold_end+1={self.scaffold_end}+1 != \
                self.next_exon.scaffold_start={self.next_exon.scaffold_start}")

    def __str__(self):
        return f"""
    Exon in gene: {self.scaffold_name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}"""
    
    def check_gene(self):
        if not self.scaffold_name == self.gene.scaffold_name:
            print(f"Intron: {self}, its gene: {self.gene}")
            raise ValueError(f"[EXON] Wrong scaffolds, gene is {self.gene.scaffold_name}, exon is {self.scaffold_name}")
        if not (self.gene.scaffold_start<=self.scaffold_start<self.scaffold_end<=self.gene.scaffold_end):
            print(f"gene {self.gene}, intron {self}")
            raise ValueError(f"""[EXON] Exon not in gene! gene={self.gene}:({self.gene.scaffold_start} {self.gene.scaffold_end}),
                            exon={self}:({self.scaffold_start} {self.scaffold_end})""")
                
    def assign_gene(self, gene):
        self.gene = gene
        self.check_gene()
    
    def merge_exons(self, next_exon_to_merge):
        print("Merging exons", self, next_exon_to_merge)
        self.sequence += next_exon_to_merge.sequence
        self.scaffold_end = next_exon_to_merge.scaffold_end
        self.next_exon = next_exon_to_merge.next_exon


class Transcript():
    def __init__(self, scaffold_name, start, end, sequence='', strand=''):
        self.scaffold_name = scaffold_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.strand = strand


def process_file(file_path):
    try:
        with open(file_path, encoding="utf-8") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                yield [int(x) if x.isnumeric() else x for x in line.split()]
    except (IOError, OSError) as exc:
        print("Error opening / processing file")
        raise exc


def create(genome_path, genes_gff_path, genes_data_type):
    '''Creates objects representing genes, introns and exons.
    Input:
        path to .fasta file,
        path to .gff/gtf file,
        data type (stringtie, gmap or manual)'''

    print("Creating genes, sequences, introns, exons and introns' class prediction")
    valid_data_type = {'stringtie', 'gmap', 'manual'}
    if genes_data_type not in valid_data_type:
        raise ValueError(f"[create()] data_type must be one of {valid_data_type}.")
    genome = read_genome(genome_path, 'fasta')
    genes = read_genes(genes_gff_path, genes_data_type)
    genes_sorted = sorted(genes.values(), key=lambda g:
        (g.scaffold_name, g.scaffold_start, g.scaffold_end))
    for gene in tqdm(genes_sorted):
        gene.add_exons(genes_data_type)
        gene.extract_sequence(genome)
        gene.create_introns()
    print("[CREATE] Genes, sequences, introns, exons created, introns' classes predicted")
    return genome, genes


# def create_manual_from_gb(genbank_path):
#     """
#     Returns genes and genome objects.
#     :param genbank_path: [
#     :return: two dictionaries, genes = {} and genome = {gene_name:scaffold_seq}
#     """
#     label_dict = create_label_dict()
#     genes = {}
#     genome = defaultdict(str)
#     with open(genbank_path) as f:
#         for record in SeqIO.parse(f, "genbank"):
#             gene_name, gene = record.name, None
#             list_of_exons = []
#             gene_start, gene_end = 1000000, 0
#             sequence = str(record.seq)
#             genome[gene_name] = sequence
#             for feature in record.features:
#                 #start, end = int(feature.location.start.position), int(feature.location.end.position)
#                 start, end = int(feature.location.start), int(feature.location.end)
#                 if feature.type == 'exon':
#                     gene_start = min(start, gene_start)
#                     gene_end = max(end, gene_end)
#                     exon = Exon(gene_name, start, end, strand='+', gene=None)
#                     list_of_exons.append(exon)
#             gene = Gene(gene_name, gene_start, gene_end, name=gene_name, strand='+', exons=[])
#             for exon in list_of_exons:
#                 exon.assign_gene(gene)
#             gene.working_exons = list_of_exons
#             gene.add_exons('manual')
#             gene.extract_sequence(genome)
#             try:
#                 gene.create_introns()
#             except Exception as exc:
#                 print("didnt create introns??? in", gene_name)
#                 print(gene)
#                 raise exc
#             for feature_listed in record.features:
#                 #start, end = int(feature_listed.location.start.position), int(feature_listed.location.end.position)
#                 start, end = int(feature_listed.location.start), int(feature_listed.location.end)
#                 if 'label' in feature_listed.qualifiers.keys():
#                     label = feature_listed.qualifiers['label'][0]
#                 elif 'standard_name' in feature_listed.qualifiers.keys():
#                     label = feature_listed.qualifiers['standard_name'][0]
#                 else:
#                     label = ''
#                 if label in label_dict:
#                     label = label_dict[label]
#                     try:
#                         intron_obj = gene.introns_dict[(start, end)]
#                     except KeyError:
#                         continue
#                     intron_obj.add_manual_annotation(label, start, end)
#             genes[gene_name] = gene
#         predict_all_introns(genes)
#     return genome, genes

def create_manual_from_gb(genbank_path, if_classifier_mixed=False):
    """
    Returns genes and genome objects.
    :param genbank_path: [
    :return: two dictionaries, genes = {} and genome = {gene_name:scaffold_seq}
    """
    label_dict = create_label_dict()
    genes = {}
    genome = defaultdict(str)
    with open(genbank_path, encoding='utf-8') as f:
        for record in SeqIO.parse(f, "genbank"):
            gene_name, gene = record.name, None
            list_of_exons = []
            gene_start, gene_end = 1000000, 0
            sequence = str(record.seq)
            genome[gene_name] = sequence
            for feature in record.features:
                start, end = int(feature.location.start), int(feature.location.end)
                if feature.type == 'exon':
                    if start < gene_start:
                        gene_start = start
                    if end > gene_end:
                        gene_end = end
                    exon = Exon(gene_name, start, end, strand='+', gene=None)
                    list_of_exons.append(exon)
            gene = Gene(gene_name, gene_start, gene_end, name=gene_name, strand='+', exons=[])
            for exon in list_of_exons:
                exon.gene = gene                
            gene.working_exons = list_of_exons
            gene.add_exons('manual')
            gene.extract_sequence(genome)
            try:
                gene.create_introns()
            except Exception as exc:
                print(gene_name)
                print(gene)
                raise exc
            for feature_listed in record.features:
                start, end = int(feature_listed.location.start), int(feature_listed.location.end)
                if 'label' in feature_listed.qualifiers.keys():
                    label = feature_listed.qualifiers['label'][0]
                elif 'standard_name' in feature_listed.qualifiers.keys():
                    label = feature_listed.qualifiers['standard_name'][0]
                else:
                    label = ''
                if label in label_dict:
                    label = label_dict[label]
                    try:
                        intron_obj = gene.introns_dict[(start, end)]
                    except KeyError:
                        continue
                    intron_obj.add_manual_annotation(label, start, end)
            genes[gene_name] = gene
        predict_all_introns(genes, if_mixed=if_classifier_mixed)
    return genes, genome

def add_manual_annotation_from_gb(gb_file, genes):
    '''To previously created genes, add manual annotations
    from a given .gb file'''
    label_dict = create_label_dict()
    count = {}
    for record in SeqIO.parse(gb_file, "genbank"):
        gene_name = record.name
        for feature_listed in record.features:
            start, end = int(feature_listed.location.start), int(feature_listed.location.end)
            if 'label' in feature_listed.qualifiers.keys():
                label = feature_listed.qualifiers['label'][0]
            elif 'standard_name' in feature_listed.qualifiers.keys():
                label = feature_listed.qualifiers['standard_name'][0]
            if label not in label_dict:
                continue
            label = label_dict[label]
            gene_obj = genes[gene_name]
            for intron_obj in gene_obj.introns:
                if (intron_obj.scaffold_start, intron_obj.scaffold_end) == (start, end):
                    intron_obj.add_manual_annotation(label, start, end)
                    if label not in count:
                        count[label] = 1
                    else:
                        count[label] += 1
                    break
    return count


def create_label_dict():
    '''Returns a dict of what manual labels translate to
    (conventional, nonconventional or intermediate intron).'''
    k_labels = ['Intron K', 'Inron K', 'Intron K (EL, EG)', 'Intron K ',
                'Exon K', "Intron K (5' partial)"]
    nk_labels = ['Intron N', 'Intron N (EH)', 'Intron NK', 'Intron NK ',
                 'Intron N ', 'Intron NK (polimorfizm)']
    i_labels = ['Intron I', 'Intron I ', 'I', 'Intron P']
    label_dict = {}
    for i in k_labels:
        label_dict[i] = 'intron_K'
    for i in nk_labels:
        label_dict[i] = 'intron_NK'
    for i in i_labels:
        label_dict[i] = 'intron_I'
    return label_dict

def create_both(genome, genes, genes_rev, data_type):
    """
    Creates genes' and genomes' object from
    two gene format files and the same genome (fasta).
    Input:
        path to genome .fasta file,
        path to first .g(f/t)f file,
        path to 2nd .g(f/t)f file,
        data type (stringtie, gmap or manual)
    """
    genome_eug, genes_eug = create(genome, genes, data_type)
    genome_eug,genes_eug_rev = create(genome, genes_rev, data_type)
    return genome_eug, genes_eug, genes_eug_rev


def read_genome(file_path, file_type):
    '''Returns a dict of records in genome.'''
    genome = defaultdict(str)
    with open(file_path, encoding='utf-8') as f:
        for record in SeqIO.parse(f, file_type):
            genome[record.id] = str(record.seq)
    return genome
    

def read_genes(filename, data_type):
    '''Redirects to a specific read_genes function
    according to data type.'''
    if data_type == 'stringtie':
        return read_genes_stringtie(filename)
    elif data_type == 'gmap':
        return read_genes_gmap(filename)
    elif data_type == 'manual':
        return read_genes_manual(filename)


def read_genes_stringtie(filename):
    genes = {}  # slownik genow
    for line in process_file(filename):
        if line[0] == '#':
            continue
        if line[2] == 'transcript':
            if line[6] in {"-", "+"}:
                gene_name = line[11].strip('";"')
                coverage = line[13].strip('";"')
                gene = Gene(line[0], line[3] - 1, line[4], name=gene_name,
                            strand=line[6], exons=[], coverage=coverage)
                genes[gene_name] = gene
        elif line[2] == 'exon':
            if not line[6] in {"-", "+"}:
                continue
            gene_name = line[11].strip('";"')
            gene_of_exon = genes[gene_name]
            exon = Exon(line[0], line[3]-1, line[4], strand=line[6], gene=gene_of_exon)
            gene_of_exon.working_exons.append(exon)
    return genes


def read_genes_gmap(filename):
    def extract_features_from_line(line):
        scf_name, start, end, strand = line[0], line[3] - 1, line[4], line[6]
        return scf_name, start, end, strand
    
    genes = {}  # slownik genow
    exons_to_add_later = []
    for line in process_file(filename):
        try:
            if line[0][0] == '#':
                continue
            if len(line) < 3:
                print(line)
                continue
            
            scf_name, start, end, strand = extract_features_from_line(line)
            if not strand in {"-", "+"}:
                continue
            
            if line[2] == 'gene':
                gene_name = search(r'ID=([\w.]+)', line[8]).groups()[0]
                gene = Gene(scf_name, start, end, name=gene_name,
                            strand=strand, exons=[])
                genes[gene_name] = gene
            elif line[2] == 'exon':
                gene_name = search(r'Parent=([\w.]+)', line[8]).groups()[0]
                gene_name = gene_name[:-5]+"path"+gene_name[-1] #parent is mrna, has to be gene
                exons_to_add_later.append((gene_name, line))
                #not adding introns here in case some occur before their mRNA (gene)

        except AttributeError as exc:
            print(line)
            raise exc

    for gene_name, line in exons_to_add_later:
        scf_name, start, end, strand = extract_features_from_line(line)
        gene_of_exon = genes[gene_name]

        exon = Exon(scf_name, start, end, strand=strand, gene=gene_of_exon)
        gene_of_exon.working_exons.append(exon)
        
    print(f"created {len(genes)} genes")
    return genes


def read_genes_manual(filename):
    genes = {}  # slownik genow
    gene = None
    genes_exons = {}
    
    for line in process_file(filename):
        if line[0] == '#':
            continue

        curr_gene_name, obj_start, obj_end, strand = line[0], line[3] - 1, line[4], line[6]
        if strand not in {"-", "+"}:
            print(line)
            return

        if line[2] == 'transcript':
            gene = Gene(curr_gene_name, obj_start, obj_end, name=curr_gene_name, strand=strand, exons=[])
            genes[curr_gene_name] = gene
            genes_exons[curr_gene_name] = []

        elif line[2] == 'exon':
            exon = Exon(curr_gene_name, obj_start, obj_end, gene=gene, strand=strand)
            genes_exons[curr_gene_name].append(exon)
            
    for gn, exs in genes_exons.items():
        g = genes[gn]
        g.working_exons = exs
        for e in exs:
            e.gene = g

    for name, g in genes.items():
        if name != g.name:
            print("co")
        for e in g.working_exons:
            if not e.gene == g:
                print(e.gene, g)
        
    return genes


def complimentary(n1, n2):
    return n1+n2 in { "AT", "CG", "GT", "TA", "GC", "TG" }


def calculate_pairing(str1, str2):
    """Returns the no. of positions where two seqs are complimentary.
    Seqs must be of the same length."""
    if len(str1) != len(str2):
        print(len(str1), len(str2))
        print(str1, str2)
        raise ValueError
    return sum( 1 for i in range(len(str1)) if complimentary(str1[i-1], str2[-i]))


def conventional_class_rate(wersja=False):
    """returns weights for introns' classes if given no. of class,
    if not returns the weights dict"""
    wagi={0:0, 1:1, 2:3, 3:5, 4:7, 5:2, 6:4, 7:6, 8:8}
    if wersja:
        return wagi.get(wersja)
    return wagi


def calculate_pyrimidine_content(seq):
    return (seq.count('C') + seq.count('T') + seq.count('Y'))/len(seq)


def compute_intron_characteristics(seq, whether_weighted_scores = True):#prev_exon_seq, intron_seq, next_exon_seq):
    '''Calculates a set of characteristics for ML predictions.
    Input can be str or Intron.'''
    if isinstance(seq, SeqIO.SeqRecord):
        seq = str(seq)
    if isinstance(seq, str):
        prev_exon_seq = seq[:5]
        intron_seq = seq[5:-5]
        next_exon_seq = seq[-5:]
    else: # type(seq)==Intron
        prev_exon_seq = seq.prev_exon.sequence[-5:]
        intron_seq = seq.sequence
        next_exon_seq = seq.next_exon.sequence[:5]
        seq = prev_exon_seq + intron_seq + next_exon_seq
    baseY = {'C', 'T'}
    baseR = {'A', 'G'}

    e_y_cnt, i_r_cnt, i_CAG_cnt, i_y_cnt, e_r_cnt = 0, 0, 0, 0, 0

    if len(prev_exon_seq)>0 and prev_exon_seq[-1] in baseY:
        e_y_cnt = 1
    #print("intron_seq: ", intron_seq)
    if len(intron_seq)<6:
        return [0, 0, 0, 0, 0, 0, 0, 0]
    if intron_seq[0] in baseR:
        i_r_cnt = 1
    if intron_seq[3:6]=="CAG" and intron_seq[-8:-5]=="CTG":
        i_CAG_cnt = 1
    if intron_seq[-1] in baseY:
        i_y_cnt = 1
    if next_exon_seq and next_exon_seq[0] in baseR:
        e_r_cnt = 1
    pl2, pl1, pl3 = intron_pairing_score(seq, whether_weighted_scores = whether_weighted_scores)
    return_array = np.array([e_y_cnt, i_r_cnt, i_CAG_cnt, i_y_cnt, e_r_cnt, pl2, pl1, pl3])
    return return_array

def intron_pairing_score(sequence, whether_weighted_scores = False):
    def pairing_length(s, start, end):
        best_pairing_length = 0.
        pairing_length = 0.
        
        #Number of pairs to iterate over
        #Needs to be limited for shorter introns
        pair_num = min(20, len(s)-start, len(s)+end)
        
        while start <= pair_num:
            if complimentary(s[start], s[end]):
                pairing_length += weigh_pairings(s[start], s[end]) if whether_weighted_scores else 1.
            else: #if the pair is not complementary
                best_pairing_length = max(best_pairing_length, pairing_length)
                pairing_length = 0.
            #przesuniecie do nastepnej pary
            start += 1
            end -= 1
        
        best_pairing_length = max(best_pairing_length, pairing_length)
        return best_pairing_length
    
    pl2 = pairing_length(sequence, 0, -3) #przesuniecie o 2
    pl1 = pairing_length(sequence, 0, -2) #przesuniecie o 1
    pl3 = pairing_length(sequence, 0, -4) #przesuniecie o 3
    return [pl2, pl1, pl3]

# # as for 17/11/22 no longer in use
# def predict_all_introns(genes):
#     print("Predicting ML classes for introns")
#     introns_to_assess, their_characteristics = [], []
#     for _, gene in list(genes.items()):
#         for intron in gene.introns:
#             introns_to_assess.append(intron)
#             their_characteristics.append(compute_intron_characteristics(intron))
#             for var in intron.variants:
#                 introns_to_assess.append(var)
#                 their_characteristics.append(compute_intron_characteristics(var))
#     print(len(introns_to_assess), len(their_characteristics))
#     print("Finished predicting ML classes for introns!")
#     assert len(their_characteristics) > 0
    
#     predictions_conv = conventional_model.predict(their_characteristics)
#     predictions_nonconv = nonconventional_model.predict(their_characteristics)
#     probas_conv = conventional_model.predict_proba(their_characteristics)
#     probas_nonconv = nonconventional_model.predict_proba(their_characteristics)
#     #@TODO tu jest cos nie tak, gdzies zwracane jest to samo dla obu model
    
#     #probas = loaded_model.decision_function(their_characteristics)
#     for (i, pred_nc, pred_c, prob_c, prob_nc) in\
#         zip(introns_to_assess, predictions_nonconv, predictions_conv, probas_conv, probas_nonconv):
#         i.ML_nonconv_score = prob_nc[1]
#         i.ML_conv_score = prob_c[0]
#         i.ML_class = (pred_c, pred_nc)
#     #for intron in introns_to_assess:
#         var_probas_conv = [v.ML_conv_score for v in i.variants]
#         var_probas_nonconv = [v.ML_nonconv_score for v in i.variants]
#         i.ML_best_conv_version = i.variants[np.argmin(var_probas_conv)] if len(var_probas_conv) else i
#         i.ML_best_nonconv_version = i.variants[np.argmax(var_probas_nonconv)] if len(var_probas_nonconv) else i
#     return

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
        # '''

def calculate_gene_ends_from_scaffold(gene_scaffold_start, gene_scaffold_end,
                               obj_scaffold_start, obj_scaffold_end, strand):
    '''
    Calculates an object's (intron/exon) start and end wrt gene start and end
    i.e. localisation of the object in the gene
    (intron.gene_start, intron.gene_end) or (exon.gene_start, exon.gene_end)
    '''
    if not gene_scaffold_start <= obj_scaffold_start < obj_scaffold_end <= gene_scaffold_end:
        raise ValueError(f'''Wrong scaffold positions
                         gene scaffold: {gene_scaffold_start}-{gene_scaffold_end})
                         obj scaffold: {obj_scaffold_start}-{obj_scaffold_end}''')

    obj_gene_start, obj_gene_end = obj_scaffold_start - gene_scaffold_start, \
        obj_scaffold_end - gene_scaffold_start
    if strand == '-':
        #if strand is reversed, gene_start and end are swapped
        obj_gene_start, obj_gene_end = obj_gene_end, obj_gene_start

    if (strand=='+' and not obj_gene_start < obj_gene_end) or \
        (strand=='-' and not obj_gene_start > obj_gene_end):
        raise ValueError(f"strand is {strand}, intron.gene_start={obj_gene_start}, intron.gene_end={obj_gene_end}")
    if not obj_scaffold_end-obj_scaffold_start == abs(obj_gene_end-obj_gene_start):
        raise ValueError(f'''Wrong sequence lengths:
                         scaffold={obj_scaffold_end}-{obj_scaffold_start}={obj_scaffold_end-obj_scaffold_start},
                         gene=abs({obj_gene_end}-{obj_gene_start})={abs(obj_gene_end-obj_gene_start)}''')
    
    return obj_gene_start, obj_gene_end

        
def predict_all_introns(genes, if_mixed=False):
    '''predicting ML classes for all introns in genome'''
    print("predict_all_introns checkpoint 1/4")
    introns_to_assess, their_characteristics = [], []
    for gene in genes.values():
        #print(gene, len(gene.introns))
        for intron in gene.introns:
            if not intron.prev_exon:
                print(f"{intron.prev_exon} || {intron} || {intron.next_exon}")
            introns_to_assess.append(intron)
            ichar = compute_intron_characteristics(intron)
            if isinstance(ichar, bool):
                print(gene, "sequence:", gene.sequence, 'ichar:', ichar)
            their_characteristics.append(ichar)
            for var in intron.variants:
                introns_to_assess.append(var)
                vchar = compute_intron_characteristics(var)
                if isinstance(vchar, bool):
                    print(var, "sequence:", var.sequence, 'vchar:', vchar)
                their_characteristics.append(vchar)
    print("predict_all_introns checkpoint 2/4")
    assert len(their_characteristics) > 0

    print("predict_all_introns checkpoint 3/4")
    if if_mixed:
        nonconv_predictions = get_nonconventional_model().predict(their_characteristics)
    else:
        predictions_conv = get_conventional_model().predict(their_characteristics)
        predictions_nonconv = get_nonconventional_model().predict(their_characteristics)
        probas_conv = get_conventional_model().predict_proba(their_characteristics)
        probas_nonconv = get_nonconventional_model().predict_proba(their_characteristics)
        
    print("predict_all_introns checkpoint 4/4")

    #probas = loaded_model.decision_function(their_characteristics)
    if if_mixed:
        for i, nc_pred in zip(introns_to_assess, nonconv_predictions):
            is_k = int(bool(i.is_conventional))
            i.ML_class = [is_k, nc_pred]
    else:
        for (i, pred_nc, pred_c, prob_c, prob_nc) in\
            zip(introns_to_assess, predictions_nonconv, predictions_conv, probas_conv, probas_nonconv):
            i.ML_nonconv_score = prob_nc[1]
            i.ML_conv_score = prob_c[0]
            i.ML_class = (pred_c, pred_nc)
    #for intron in introns_to_assess:
        var_probas_conv = [v.ML_conv_score for v in i.variants]
        var_probas_nonconv = [v.ML_nonconv_score for v in i.variants]
        i.ML_best_conv_version = i.variants[np.argmin(var_probas_conv)] if len(var_probas_conv) else i
        i.ML_best_nonconv_version = i.variants[np.argmax(var_probas_nonconv)] if len(var_probas_nonconv) else i

def weigh_pairings(nn1: str, nn2: str) -> float:
    return {
            "AT": 0.5,
            "TA": 0.5,
            "GC": 1.0,
            "CG": 1.0,
            "GT": 0.375,
            "TG": 0.375,
        }.get(nn1+nn2, 0.0)


def get_introns_demulti(genes, if_halfway=False):
    """
    Get lists of conventional and nonconventional introns without repetitions
    from manual annotations.

    Args:
        genes (dict): _description_
        if_halfway (bool): returns lists of K, NK, I introns @TODO figure out why it is used

    Returns:
        introns_K_demulti (list): a list of unique conventional introns
        introns_NK_demulti (list): a list of unique nonconventional introns
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
    
    if if_halfway:
        return introns_NK, introns_K, introns_I

    def demultiply_group(introns_X):
        introns_X_demulti = []
        single, double = [], {}
        for intron in introns_X:
            if intron.sequence in single:
                for i, single_i in enumerate(single):
                    if intron.sequence != single_i:
                        continue
                    name_1, name_2 = str(introns_X_demulti[i]), str(intron)
                    pos_1, pos_2 = name_1.split()[1:], name_2.split()[1:]
                    if pos_1 == pos_2:
                        continue
                    if name_1 not in double:
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


def liczenie_intronow(name, genome_eug, genes_eug, file_type=None, czy_fasta_gtf=False, czy_wykresy=False):
    if czy_fasta_gtf:
        genome,genes=create(genome_eug, genes_eug, file_type)
        genome_eug=genome
        genes_eug=genes
    conv_count=0
    nconv_count=0
    both_count=0
    all_count=0
    non_count=0
    conventional_classes=[0,0,0,0,0,0,0,0]
    nonconventional_classes=[0,0,0,0,0,0,0,0,0]

    for _, gene in list(genes_eug.items()):
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
    print(f"""\nKonwencjonalne: {conv_count}
          niekonwencjonalne: {nconv_count}
          oba: {both_count}
          inne: {non_count}
          wszystkie: {all_count}""")
    stats=[name, conventional_classes, nonconventional_classes,
           conv_count,nconv_count, both_count,non_count,all_count]
    if czy_wykresy:
        wykresy_liczenia(stats)
    return stats


def wykresy_liczenia(stats, wykresy=("K"), tabele=True, save_wykresy=False):
    if wykresy:
        colors=[np.random.rand(3,) for i in range(len(stats))]
        for i in wykresy:
            if i not in ["K", "NK", "KNBN"]:
                print("Wrong plot values")
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
        for i, d in enumerate(data):
            ax.bar(X+width*i, d, color=colors[i], width=width)
        ax.set_title('Klasy intronow konwencjonalnych')
        ax.set_xlabel('Klasy konwencjonalne')
        ax.set_ylabel('% intronów')
        ax.set_xticks(X)
        ax.set_xticklabels(('Class 1', 'Class 2', 'Class 3', 'Class 4',
                            'Class 5', 'Class 6', 'Class 7', 'Class 8'))
        ax.legend(legenda)
        if save_wykresy:
            plt.savefig("klasy_konwencjonalne.png")
    
    if "NK" in wykresy:
        lista1=[i[2] for i in klasy]
        data=[[j*100 for j in l] for l in lista1]
        legenda=[i[0] for i in klasy]
        Y = np.arange(9)
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_axes([0,0,1,1])
        width = 1/(len(data)+1)
        for i, d in enumerate(data):
            ax.bar(Y+width*i, d, color=colors[i], width=width)
        ax.set_title('Klasy intronów niekonwencjonalnych')
        ax.set_xlabel('Klasy niekonwencjonalne')
        ax.set_ylabel('% intronów')
        ax.set_xticks(Y)
        ax.set_xticklabels(('Class 1', 'Class 2', 'Class 3', 'Class 4', 'Class 5', 'Class 6', 'Class 7', 'Class 8', 'Class 9'))
        ax.legend(legenda)
        if save_wykresy:
            plt.savefig("klasy_niekonwencjonalne.png")
            
        if "KNBN" in wykresy:
            data=[[j*100 for j in l[1:-1]] for l in counts]
            #print(data)
            labels=[i[0] for i in counts]
            Z = np.arange(4)
            fig = plt.figure(figsize=(7,5))
            ax = fig.add_axes([0,0,1,1])
            width = 1/(len(data)+1)
            for i, d in enumerate(data):
                ax.bar(Z+width*i, d, color=colors[i], width=width)
            ax.set_title('Ilości wszystkich intronów')
            ax.set_xlabel('Przynależność do klas')
            ax.set_ylabel('% intronów')
            ax.set_xticks(Z)
            ax.set_xticklabels(('Konwencjonalne','Niekonwencjonalne', 'Oba','Żadne'))
            ax.legend(labels)
        if tabele:
            print(tabulate(stats, headers=['nazwa', 'klasy konw.', 'klasy niekonw.', 'konwencjonalne', 'niekonwencjonalne', 'oba', 'żadne', "wszystkie"]))
        return
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


####    Model management

def load_model(fn: str):
    """
    Load a pickled sklearn model from filename
    """
    with open(fn, 'rb') as fd:
        out = pickle.load(fd)
    return out


#Default models to use
models = {
    'N': load_model("./Models/29_11_K_model.sav"),
    'C': load_model("./Models/29_11_NK_model.sav")
}


def get_conventional_model():
    '''Gives the currently used model for conventional introns.'''
    return models['C']

def get_nonconventional_model():
    '''Gives the currently used model for nonconventional introns.'''
    return models['N']

