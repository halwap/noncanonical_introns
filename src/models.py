import pickle

nonconventional_model = pickle.load(open('../../Models/29_11_K_model.sav', 'rb'))
conventional_model = pickle.load(open('../../Models/29_11_NK_model.sav', 'rb'))
#nonconventional_model = pickle.load(open('./files_for_classifiers/29_11_NK_model.sav', 'rb'))
#conventional_model = pickle.load(open('./files_for_classifiers/29_11_K_model.sav', 'rb'))
#nonconventional_model = pickle.load(open('./files_for_classifiers/15_11_NK_model.sav', 'rb'))
#conventional_model = pickle.load(open('./files_for_classifiers/15_11_K_model.sav', 'rb'))


def get_conventional_model():
    '''Gives the currently used model for conventional introns.'''
    return conventional_model

def get_nonconventional_model():
    '''Gives the currently used model for nonconventional introns.'''
    return nonconventional_model