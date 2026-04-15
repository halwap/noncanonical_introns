import os
import sys
import src.introns as introns
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import csv
from tabulate import tabulate
import subprocess


wykresy = ["K", "NK", "KNBN"]
czy_tabele=True
save_wykresy = False
args=sys.argv
tekstdoliczenia='TekstDoLiczenia_maly.txt'
struktury=False

for arg in range(1, len(args)):
    if args[arg]=="-wykresy":
        wykresy=[]
        i=arg
        while len(args)>i+1 and args[i+1] in ["K", "NK", "KNBN"]:
            wykresy.append(args[i+1])
            i+=1
    elif args[arg]=='-tabele':
        czy_tabele=True
    elif args[arg]=="-savewykresy":
        save_wykresy=True
    elif args[arg]=='-plik':
        tekstdoliczenia=args[arg+1]
    elif args[arg]=='-struktury':
        struktury=True
print("wykresy: %s, tabele: %s, zapisać tabele: %s, plik do licznenia: %s, struktury: %s" %(" ".join(wykresy), czy_tabele, save_wykresy, tekstdoliczenia, struktury))

pliki=[l.strip().split(", ") for l in open(tekstdoliczenia, 'r').readlines()]

ggs=[]
s=[]
for name, genome_eug, genes_eug in pliki:
    print(name)
    genome, genes = introns.create_genes(genome_eug, genes_eug)
    ggs.append((name, genome,genes))
    stats = introns.liczenie_intronow(name, genome_eug, genes_eug, czy_fasta_gtf = True)
    s.append(stats)
    print("creating %s done" %name, list(genes.items())[0][1].introns[0].is_struct_nonconv_1)
#print(s)
exit()

introns.wykresy_liczenia(s, wykresy = wykresy, tabele = czy_tabele, save_wykresy = save_wykresy)
if not struktury:
    exit()


    
def wszystkie_introny_wersje_fasta(name, genome, genes, path='', dlugosckoncow=25):
    #genome,genes=introns.create_genes(genome_eug, genes_eug)
    #print("creating done")    
    nazwapliku=path+'%s_25zmarginesami.fasta' %name
    print(nazwapliku)
    with open(nazwapliku, 'w') as f:
        f.write(";nazwa - nazwagenu;scaffold start - scaffold end;wersja(podst=0);klasa konw;klasa niekonw\n")
        f.write(";seq - exon[-5:] intron[:%d] intron[-%d:] exon[:5]\n" %(dlugosckoncow, dlugosckoncow))
        print("headline written")
        for nazwa, gene in list(genes.items()):
            #print(nazwa)
            for intron in gene.introns:
                s=intron.sequence
                vconv = intron.is_conventional if intron.is_conventional else 0
                vnonconv = intron.is_nonconventional
                #vstruct = 1 if intron.is_struct_nonconv_1 else 2 if intron.is_struct_nonconv_2 else 0
                f.write(">%s;%s-%s;0;%s;%s\n" %(nazwa,intron.scaffold_start, intron.scaffold_end, vconv,vnonconv))
                f.write(intron.prev_exon.sequence[-5:]+s[:dlugosckoncow]+"AAAAAAAAAA"+s[-dlugosckoncow:]+intron.next_exon.sequence[:5]+"\n")

                for index, var in enumerate(intron.variations):
                    s=var.sequence
                    vconv = var.is_conventional if intron.is_conventional else 0
                    vnonconv = var.is_nonconventional
                    f.write(">%s;%s-%s;%d;%s;%s\n" %(nazwa, var.scaffold_start, var.scaffold_end, index,vconv,vnonconv))
                    f.write(var.prev_exon.sequence[-5:]+s[:dlugosckoncow]+"AAAAAAAAAA"+s[-dlugosckoncow:]+var.next_exon.sequence[:5]+"\n")
    return nazwapliku

def count_nconv_kilkawersji(listaplikow, definicja_nonconv="seq", path=""):
    ulamki=[]
    dane_liczbowe=[]
    for plik in listaplikow:
        name='_'.join(plik.split('_')[:2])    
        #if plik.startswith("wszystkie_introny_wersje_E."):
        #    name=plik.split("_")[3]
        #elif plik.startswith("generowane_introny_"):
        #    name=plik.split(".")[0]
        #else:
        #    print("zla nazwa pliku")
        #    break
        
        plik=".".join(plik.split(".")[:-1])
        plikfasta, plikdbn=plik+".fasta", plik+".dbn"
        fasta=open(plikfasta, 'r').readlines()
        if fasta[0][0]==";": fasta=fasta[2:]
        dbn = open(plikdbn, 'r').readlines()
        
        fasta = [ x+y for x,y in zip(fasta[0::2], fasta[1::2]) ]
        dbn = [x+y+z for x,y,z in zip(dbn[0::3], dbn[1::3], dbn[2::3])]
        
        l=len(fasta) if len(fasta)==len(dbn) else None
        
        global conv_count, nconv_count, both_count, all_count, non_count
        conv_count, nconv_count, both_count, all_count, non_count = 0,0,0,0,0
        prev_var = [0, False, False] #wersja, konwencjonalny, niekonwencjonalny struktur

        def dodaj(prev_var):
            global conv_count, nconv_count, both_count, all_count, non_count
            var, conv_ver, nonconv_str=prev_var
            all_count+=1

            if nonconv_str and conv_ver: both_count+=1
            elif nonconv_str and not conv_ver: nconv_count+=1
            elif not nonconv_str and conv_ver: conv_count+=1
            else: non_count+=1
            
        for i in range(l):
            f, d = fasta[i].split("\n"), dbn[i].split("\n")
            if f[1] != d[1].replace("U", "T"):
                print("zle sekwencje")
                break
            
            var,conv_ver = f[0].split(";")[2:4]
            var=int(var)
            conv_ver=int(conv_ver)
            vienna = d[2].split(" ")[0]
            seq=d[1]
            
            if definicja_nonconv=="struct":
                #warunki na niekonwencjonalne
                nonconv_str=0
                nazwa_pliku_png = "test_parowaniakons_przesuniecieo2_Ypoczkon"
                if not '(' in vienna[-25:-5] and not ')' in vienna[5:25]:
                    #if vienna.count('(', 0, 25)>5:
                    #    nonconv_str = 2
                    #    if seq[-5] in ["G", "A"]:
                    #        nonconv_str = 1
                    parowania_kons = vienna[8:11]=="(((" and vienna[-13:-10]==")))" #pozycje 4,5,6 (-8,-7,-6) -> indeksy 3,4,5 (-8,-7,-6) -> z 5 nt eksonu 8,9,10 (-13,-12,-11)
                    parowan_5 = vienna.count('(', 0, 25)>=5
                    przesuniecie = vienna[:8].count("(") == vienna[-10:].count(")")
                    puryna_pocz = seq[5] in ["G", "A"]
                    puryna_koniec =  seq[-5] in ["G", "A"]
                    pirymidyna_pocz = seq[4] in ["U", "C", "T"]
                    pirymidyna_koniec = seq[-4] in ["U", "C", "T"]
                    
                    zestaw_zasad = [parowania_kons, przesuniecie, pirymidyna_koniec, pirymidyna_pocz]
                    if False not in zestaw_zasad:
                        nonconv_str = 1
                    #nazwa_pliku_png = ''
                #warunki na niekonwencjonalne

            if var and var>prev_var[0]:
                if conv_ver>0: prev_var[1] = True
                if nonconv_str>0: prev_var[2] = True
                prev_var[0]=var
            elif not var:
                dodaj(prev_var)
                prev_var=[0,bool(conv_ver),bool(nonconv_str)]
        ulamki.append([name, nconv_count/all_count, all_count])
        dane_liczbowe.append([name, nconv_count, all_count])
    for count in ulamki:
        print(count)
    for i in dane_liczbowe:
        print(i)        
    return ulamki, dane_liczbowe#, nazwa_pliku_png

path_working="./ViennaRNA-2.4.18/working_fastas/"
#sztuczne=["generowane_introny_0.fasta", "generowane_introny_1.fasta", "generowane_introny_2.fasta", "generowane_introny_3.fasta", "generowane_introny_4.fasta"]
gat=[]
for name, genome,genes in ggs:
    nazwapliku=wszystkie_introny_wersje_fasta(name, genome,genes,path=path_working)
    gat.append(nazwapliku)
print("zrobione fasty")

for plik in gat:
    p=plik.split('/')[-1].split('.')[0]
    command = ['/usr/local/bin/RNAfold', '--noPS',  '--command=constr.txt', '-j8', '-i', path_working+p+'.fasta']
    #print(command)
    #print(p+'.dbn')
    with open(path_working+p+'.dbn', 'w') as output_file:
        _p = subprocess.run(command, stdout=output_file)
    #print(plik)

definiowanie="struct"
#path="/home/semik/projekty/noncanonical_introns/fasty_i_dbn/"
#ulamki, dane_liczbowe, nazwa_pliku_png = count_nconv_kilkawersji(gat, definiowanie, path=path_working)
ulamki, dane_liczbowe = count_nconv_kilkawersji(gat, definiowanie, path=path_working)
print("liczby", tabulate(dane_liczbowe))
exit()
mean_list=mean_generowanych(ulamki)

colors = ['tab:blue', 'tab:green', 'tab:red', 'darkorchid']
podpis="a+b+c. Min. 5 parowań, w tym na konserwowanych pozycjach z przesunięciem o 2"
#wykres_licznosci_maly(ulamki[:3]+[mean_list], colors, podpis, definiowanie=definiowanie, path=path_working, nazwa_pliku_png = nazwa_pliku_png)

#/home/semik/projekty/noncanonical_introns/ViennaRNA-2.4.18/working_fastas/bugtest_25zmarginesami.fasta 
#/home/semik/projekty/noncanonical_introns/ViennaRNA-2.4.18/working_fastas
#/home/semik/projekty/noncanonical_introns/ViennaRNA-2.4.18/working_fastas/bugtest_25zmarginesami.fasta 


'''def wykres_licznosci_maly(counts, colors, podpis, definiowanie, path="", nazwa_pliku_png="liczby_intronow"):
    nazwa_pliku_png = nazwa_pliku_png + "_%s.png" %definiowanie
    print(counts)
    data=[l[1]*100  for l in counts]
    labels=[i[0] for i in counts]
    print(labels, data)
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0,0,1,1])
    
    ax.bar(x=labels, height=data,color=colors)
    ax.set_xlabel(podpis)
    ax.set_ylabel('% intronów')
    plt.ylim([0,12.5])
    plt.show()
    #fig.savefig(path+nazwa_pliku_png , bbox_inches='tight')'''
