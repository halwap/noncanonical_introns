import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer


def cut_junctions(cut, sequences):
    '''Preprocess data in different ways by cutting off the conventional or all splices.
    0 - dont cut, 1 - only conventional, 2 - both'''
    new_seqs = []
    for i in sequences:
        if cut == 0:
            new_seqs.append(i)
            continue
        left_anchor = i[1][3:5]
        right_anchor = i[1][-5:-3]
        if left_anchor == 'GT' and right_anchor == 'AG' or left_anchor == 'CT' and right_anchor == 'AC':
            new_seqs.append(i[1][5:-5])
        else:
            if cut == 1:
                new_seqs.append(i)
            else:
                new_seqs.append(i[1][8:-9])
    return new_seqs


def kmers(sequence, k):
    return [sequence[i: i + k] for i in range(len(sequence)-1)]


def split(k, seqs):
    return [" ".join(kmers(sequence, k)) for sequence in seqs]


def preprocess(cut, k, sequences, file_out):
    """Optionally cut off junctions, then split into k-length oligonucleotides and convert into TF_IDF."""
    sequences_cut = cut_junctions(cut, [x[1] for x in sequences])
    split_seqs = split(k, sequences_cut)
    tfidf_vectorizer = TfidfVectorizer(use_idf=True, lowercase=False)
    tfidf_vectorizer_vectors = tfidf_vectorizer.fit_transform(split_seqs)
    print(type(tfidf_vectorizer_vectors))
    np.save(file_out, tfidf_vectorizer_vectors)
