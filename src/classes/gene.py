from bisect import insort
from copy import copy
import numpy as np
from src.tools import reverse_complement, calculate_gene_ends_from_scaffold, compute_intron_characteristics
from src.models import nonconventional_model, conventional_model
from .genomic_sequence import GenomicSequence
from .transcript import Transcript
from .intron import Intron

class Gene(GenomicSequence):
    '''
    Gene class representing a gene with exons and introns.
    '''
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

    def __str__(self):
        return f"""
    Gene {self.name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}"""
    
    def append_introns(self, intron):
        '''Adds the given intron to the list of the gene's introns'''
        self.introns.append(intron)
    
    # def create_intron_structures(self, path_working="./ViennaRNA-2.4.18/working_fastas/", path_RNAfold = '/usr/local/bin/'):
    #     '''To be used with RNAfold'''
    #     p=plik.split('/')[-1].split('.')[0]
    #     command = [path_RNAfold+'RNAfold', '--noPS',  '--command=constr.txt', '-j8', '-i', path_working+p+'.fasta']
    #     with open(path_working+p+'.dbn', 'w') as output_file:
    #         _p = subprocess.run(command, stdout=output_file)

    def append_exons(self, exon):
        '''Adds exon to the exons list in the correct order'''
        insort(self.exons, exon, key=lambda _exon: _exon.scaffold_start)
        self.exons.append(exon)


    def add_exons(self, genes_data_type, print_info=False):
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
            if print_info:
                print(f"[loop start] prev_exon {exon_s_e(prev_exon)} curr_exon {exon_s_e(curr_exon)}")
            if prev_exon:
                if genes_data_type == 'manual' and \
                    prev_exon.scaffold_end == curr_exon.scaffold_start:
                    #merging exons if there is no intron between them
                    new_exon = copy(prev_exon)
                    new_exon.sequence = prev_exon.sequence + curr_exon.sequence
                    new_exon.scaffold_end = curr_exon.scaffold_end
                    new_exon.next_exon = curr_exon.next_exon
                    if print_info:
                        print(f"merging {exon_s_e(prev_exon)} with {exon_s_e(curr_exon)} into {exon_s_e(new_exon)}")
                    if len(self.exons)>0:
                        self.exons[-1] = new_exon
                    else:
                        self.exons.append(new_exon)
                    prev_exon, curr_exon = new_exon, new_exon.next_exon
                else:
                    prev_exon.next_exon = curr_exon
                    curr_exon.prev_exon = prev_exon
                    self.exons.append(curr_exon)
                    if print_info:
                        print(f"appended {exon_s_e(curr_exon)}")
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
            self.sequence = reverse_complement(sequence)
            self.expanded_sequence = reverse_complement(expanded_sequence)
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
            nonconv_predictions = nonconventional_model.predict(their_characteristics)
            for i, nc_pred in zip(introns_to_assess, nonconv_predictions):
                is_k = int(bool(i.is_conventional))
                i.ML_class = [is_k, nc_pred]
        else:
            # WAZNE nie pomylic indeksow, conv to 1 a nonconv to 0
            conv_predictions = 1-conventional_model.predict(their_characteristics)
            conv_scores = conventional_model.predict_proba(their_characteristics)[:, 1]
            nonconv_scores = nonconventional_model.predict_proba(their_characteristics)[:, 0]
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
