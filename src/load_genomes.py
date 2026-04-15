
from collections import defaultdict
from re import search
from Bio import SeqIO
from tqdm import tqdm
import numpy as np

from src.classes.exon import Exon
from src.classes.gene import Gene

from src.tools import compute_intron_characteristics
from src.models import nonconventional_model, conventional_model


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
    if data_type == 'gmap':
        return read_genes_gmap(filename)
    if data_type == 'manual':
        return read_genes_manual(filename)
    print("Wrong data type:", data_type)


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



def predict_all_introns(genes, if_mixed=False):
    '''predicting ML classes for all introns in genome'''
    print("predict_all_introns checkpoint 1/4")
    introns_to_assess, their_characteristics = [], []
    for _, gene in list(genes.items()):
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
        nonconv_predictions = nonconventional_model.predict(their_characteristics)
    else:
        predictions_conv = conventional_model.predict(their_characteristics)
        predictions_nonconv = nonconventional_model.predict(their_characteristics)
        probas_conv = conventional_model.predict_proba(their_characteristics)
        probas_nonconv = nonconventional_model.predict_proba(their_characteristics)
        
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
    for i in introns_to_assess:
        var_probas_conv = [v.ML_conv_score for v in i.variants]
        var_probas_nonconv = [v.ML_nonconv_score for v in i.variants]
        i.ML_best_conv_version = i.variants[np.argmin(var_probas_conv)] if len(var_probas_conv) else i
        i.ML_best_nonconv_version = i.variants[np.argmax(var_probas_nonconv)] if len(var_probas_nonconv) else i



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

def HELP_load_default_genome_genes(species=None, with_reversed=False, if_classifier_mixed=False):
    '''
    helping function, not to be used eventually
    '''
    path_a = "../../domowe_introny/dane_wejsciowe/"
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
        return

    if not with_reversed:
        r = create(genome_file, genes_file, "stringtie")
    else:
        r = create_both(genome_file, genes_file, genes_rev_file, "stringtie")

    return r #genome_eug, genes_eug, [genes_eug_rev]
