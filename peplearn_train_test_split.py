from pandas import Series
from utils import load_data
from iedb import mhc_datasets

def get_matching_data(mhc, convolved):
    '''
    Given a dataset of peptides and a dataset of peptides convolved
    onto alleles, select the peptides in the second dataset which are
    present in the first.
    '''
    indices = []
    for index, allele in enumerate(mhc[2]):
        if allele.split('-')[0] == 'HLA':
            indices.append(index)
    seqs = []
    for index in indices:
        seq = (convolved['peptide'] == mhc[0][index])
        allele = (convolved['mhc'] == mhc[2][index])
        seqs += list(convolved[seq & allele].index)
    return Series(seqs)

def generate_tts_datasets(data_loc, arr_cols=None, names=None):
    '''
    Given the location of a dataset splits the dataset into train
    test and split files based on whether each peptide is train test or split
    in the peplearn dataset. Saves the split dataset to csvs.

    Arguments -------
    data_loc: location of dataset to split
    arr_cols: optional, name of columns of dataset with numpy arrays as strings
    names: optional, names to give converted columns
    '''
    tables = ['mhc_test1', 'mhc_test2', 'mhc_train']
    convolved_mhc = load_data(data_loc, arr_cols=arr_cols, names=names)
    datasets = []
    for table in tables:
        peps = get_matching_data(mhc_datasets(table), convolved_mhc)
        dataset = convolved_mhc.loc[peps].copy()
        dataset.reset_index(drop=True).to_csv(table + '_' + data_loc + '.csv')