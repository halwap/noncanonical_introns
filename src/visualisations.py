import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate

from src.load_genomes import create
from src.tools import intron_pairing_score

warnings.filterwarnings("ignore",
                        #message="divide by zero encountered in divide",
                        category=UserWarning)


def counting_introns(name, genome_eug, genes_eug, file_type=None, czy_fasta_gtf=False, czy_wykresy=False):
    '''
    Counts introns from genome+genes:
        All - how many there are total,
        Total of conventional, nonconvential introns
        No. of introns of each class in conv, nonconv
        Both - introns that classify both as conv and nonconv (intermediate?)
        None - introns that were classified neister as conv nor nonconv
    Optionally prints the values.
    Optionally plots the values.
    Returns a list consisting of:
        [Name, [no. of each class in conventional],
               [no. of each class in nonconventional],
               All conventional, all Nonconventional,
               Intermediate (both conv and nonconv),
               Other intrins,
               No. of all introns]
    '''
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
    print(f"""\t{name}
          Conventional: {conv_count}
          Nonconventional: {nconv_count}
          Both: {both_count}
          Other: {non_count}
          All: {all_count}""")
    stats=[name, conventional_classes, nonconventional_classes,
           conv_count,nconv_count, both_count,non_count,all_count]
    if czy_wykresy:
        plot_class_counts(stats)
    return stats

def plot_class_counts(stats, wykresy=("K"), tabele=True, save_wykresy=False):
    colors=[]
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

# not used for now
# def create_genes(genome,genes):
#     genome_eug = read_genome(genome)
#     genes_eug = read_genes(genes)
#     for name, gene in list(genes_eug.items()):
#         gene.extract_sequence(genome_eug)
#         gene.create_introns()
#     return genome_eug,genes_eug


def make_heatmap_of_pairings(introns_set, nk_or_k, which_pls_tps="tp", shifts=(2, 1), score_weighed=False):
    '''
    introns_set - preferably is the output of get_introns_demulti
    which_pls_tps:
        tp - total pairings
        pl - pairing lengths
        [tp2, tp1] - total pairings in 2bp shift vs 1bp shift
    '''
    
    assert nk_or_k in ["k", "nk"]
    assert isinstance(score_weighed, bool)
    assert which_pls_tps in ["tp", "pl"]
    
    #print(f"drawing for {nk_or_k}, score_weighed={score_weighed}, {which_pls_tps}")
    score_columns=[f"{which_pls_tps}{s}" for s in shifts]

    def get_introns_list_scores(sequences, whether_weighted_scores = False):
        return pd.DataFrame(np.array([([i]+intron_pairing_score(i.sequence, whether_weighted_scores = whether_weighted_scores)) \
            for i in sequences]), columns=["intron_object", "pl2", "pl1", "pl3", "tp2", "tp1", "tp3"])

    if nk_or_k=="nk":
        introns_scores_set = get_introns_list_scores(introns_set, score_weighed)
    else:
        introns_scores_set = get_introns_list_scores(introns_set, score_weighed)

    pls_tps = introns_scores_set[score_columns]
    pls_tps_df = pd.DataFrame(0., index = np.arange(21), columns = np.arange(21))
    # df z zerami, za każdą parę (x=parowanie w 1/3, y=parowanie w 2) dodajemy 1 do df'a
    for _, (x, y) in pls_tps.iterrows():
        #print(x,y)
        pls_tps_df.iloc[int(x),int(y)] += 1

    sns.heatmap(pls_tps_df, annot=True)
    set_name = "Nonconventional" if nk_or_k=="nk" else "Conventional"
    score_name = "total pairing #" if which_pls_tps=="tp" else "pairing lengths"
    if score_weighed:
        score_name = "weighed "+score_name
    else:
        score_name = "non-weighed "+score_name

    plt.title(f'{set_name} {score_name}: {shifts[0]} vs {shifts[1]}')
    plt.xlabel(f"shift by {shifts[1]}")
    plt.ylabel(f"shift by {shifts[0]}")
    plt.show()

def draw_all_heatmaps_of_pairings(introns_nk_set=None, introns_k_set=None, score_weighed=None, score_type=None):
    #rekurencyjne
    assert score_type in ["tp", "pl", None]
    assert introns_nk_set or introns_k_set

    score_weighed = [score_weighed] if score_weighed is not None else [True, False]
    score_type = [score_type] if score_type else ["tp", "pl"]
    shifts = [(2,1), (2,3)]
    
    #print("score_weighed, score_type, shifts:", score_weighed, score_type, shifts)
    
    if introns_nk_set:
        for s in shifts:
            for s_w in score_weighed:
                for s_t in score_type:
                    make_heatmap_of_pairings(introns_nk_set, "nk", which_pls_tps=s_t, shifts=s, score_weighed=s_w)
                
    if introns_k_set:
        for s in shifts:
            for s_w in score_weighed:
                for s_t in score_type:
                    make_heatmap_of_pairings(introns_k_set, "k", which_pls_tps=s_t, shifts=s, score_weighed=s_w)

def add_value_labels(xs,ys,ax):
    for x,y in zip(xs,ys):
        ax.text(x,y,y)
    
def plot_correct_positions(introns_set, title=None):
    fig, ax = plt.subplots(1, 2, figsize=(17,5))
    fig.suptitle(title)
    labels = ["1st pos.\ncorrect\n(classified\nas xK)",
              "other positions\ncorrect\n(not classified\nas xK)",
              "only 1st\nclassified\nas xK",
              "no positions\nclassified\nas xK",
              "first pos.\ncorrect,\nother wrong",
              "if any\nposition\nis not xK",
              "all"]
    nk_pos=[0,0,0,0,0,0,0]
    k_pos=[0,0,0,0,0,0,0]
    
    for i in introns_set:
        if i.man_annotation=="intron_NK":
            v1,v2=0,0
            if i.ML_class[1]:
                v1=1
                nk_pos[0]+=1
            if 1 not in [v.ML_class[1] for v in i.variants]:
                v2=1
                nk_pos[1]+=1
            if v1 and v2:
                nk_pos[2]+=1
            if not v1 and not v2:
                nk_pos[3]+=1
            if v1 and not v2:
                nk_pos[4]+=1
            if 1 in [v.ML_class[0] for v in i.variants]:
                nk_pos[5]+=1
            nk_pos[-1]+=1
        elif i.man_annotation=="intron_K":
            v1,v2=0,0
            if i.ML_class[0]:
                v1=1
                k_pos[0]+=1
            if 1 not in [v.ML_class[0] for v in i.variants]:
                v2=1
                k_pos[1]+=1
            if v1 and v2:
                k_pos[2]+=1
            if not v1 and not v2:
                k_pos[3]+=1
            if v1 and not v2:
                k_pos[4]+=1
            if 1 in [v.ML_class[1] for v in i.variants]:
                k_pos[5]+=1
            k_pos[-1]+=1
        else:
            #print(i.man_annotation)
            pass
    
    ax[0].set_title("Introns manually annotated as NK")
    ax[0].bar(labels, nk_pos)
    ax[1].set_title("Introns manually annotated as K")
    ax[1].bar(labels, k_pos)
    add_value_labels(labels, nk_pos, ax[0])
    add_value_labels(labels, k_pos, ax[1])
    plt.show()

def plot_man_annotation_positions_counts(introns_set, title_add=""):
    plt.figure(figsize=(4,3))
    plt.title("No. of positions classified as noncanonical "+title_add)
    nk_pos_in_nk=[]
    nk_pos_in_k=[]
    for i in introns_set:
        l = [v.ML_class[1] for v in i.variants]+[i.ML_class[1]]
        if i.man_annotation == "intron_NK":
            nk_pos_in_nk.append(sum(l))
        elif i.man_annotation == "intron_K":
            nk_pos_in_k.append(sum(l))
    bins = list(range(0, max(nk_pos_in_k+nk_pos_in_nk)+1))
    plt.hist([nk_pos_in_k, nk_pos_in_nk], bins=bins, label=["in introns manually\nannotated as K", "in introns manually\nannotated as NK"])
    plt.xticks(bins[:-1])
    plt.legend()
    plt.show()
