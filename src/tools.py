from Bio import SeqIO
import numpy as np

def complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    letters = [complement_dict[base] for base in seq]
    return ''.join(letters)


def reverse_complement(seq):
    return complement(seq[::-1])

def complimentary(n1, n2):
    return {n1, n2} in [{'A', 'T'}, {'C', 'G'}, {'G', 'T'}]

def calculate_pairing(str1, str2):
    """Returns the no. of positions where two seqs are complimentary.
    Seqs must be of the same length."""
    if len(str1) != len(str2):
        print(len(str1), len(str2))
        print(str1, str2)
        raise ValueError
    counter = 0
    for i in range(len(str1)):
        if complimentary(str1[i-1], str2[-i]): counter += 1
    return counter


def conventional_class_rate(wersja=False):
    """returns weights for introns' classes if given no. of class,
    if not returns the weights dict"""
    wagi={0:0, 1:1, 2:3, 3:5, 4:7, 5:2, 6:4, 7:6, 8:8}
    if wersja:
        return wagi.get(wersja)
    return wagi

def calculate_pyrimidine_content(seq):
    count = 0
    for char in seq:
        if char in 'CTY':
            count += 1
    return count/len(seq)

def weigh_pairings(nn1, nn2):
    nn = nn1+nn2
    if nn in ["AT", "TA"]: return 0.5
    if nn in ["GC", "CG"]: return 1
    if nn in ["GT", "TG"]: return 0.375
    return 0

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

def compute_intron_characteristics(seq, whether_weighted_scores = True):#prev_exon_seq, intron_seq, next_exon_seq):
    '''Calculates a set of characteristics for ML predictions.
    Input can be str or Intron.
    
    Checked characteristics are:
    [   0: eY | i
        1: e | Ri
        2: e | nnn CAG ... CTG nnnnn | e
        3: iY | e
        4: i | Re
        5-7: pairing score
    ]
    '''

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
    pl2, pl1, pl3 = intron_pairing_score(seq, whether_weighted_scores = whether_weighted_scores)[:3]
    return_array = np.array([e_y_cnt, i_r_cnt, i_CAG_cnt, i_y_cnt, e_r_cnt, pl2, pl1, pl3])
    return return_array

def intron_pairing_score(sequence, whether_weighted_scores = False):
    '''
    @TODO sprawdzic czy ta funkcja robi co ma robic, bo wyglada podejrzanie
    '''
    def pairing_length(s, start=0, end=-1):
        n = len(s[start:end])/2
        x, y = start, end
        best_pairing_length = 0.
        pairing_length = 0.
        total_pairings = 0.
        while x<=20:
            if total_pairings > 20:
                print("x =", x, total_pairings)
            if complimentary(s[x], s[y]):
                #print(s[x], s[y], introns.weigh_pairings(s[x], s[y]))
                pairing_length += 1. if whether_weighted_scores is False else weigh_pairings(s[x], s[y])
                total_pairings += 1. if whether_weighted_scores is False else weigh_pairings(s[x], s[y])
            else: #if the pair is not complementary
                best_pairing_length = max(best_pairing_length, pairing_length)
                pairing_length = 0.
            #przesuniecie do nastepnej pary
            x += 1
            y -= 1
        best_pairing_length = max(best_pairing_length, pairing_length)
        # dlaczego to tu jest? co to mialo robic??
        # if whether_weighted_scores is True:
        #     best_pairing_length, total_pairings = best_pairing_length, total_pairings
        return best_pairing_length, total_pairings

    pl2, tp2 = pairing_length(sequence, 0, -3) #przesuniecie o 2
    pl1, tp1 = pairing_length(sequence, 0, -2) #przesuniecie o 1
    pl3, tp3 = pairing_length(sequence, 0, -4) #przesuniecie o 3
    return [pl2, pl1, pl3, tp2, tp1, tp3]

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