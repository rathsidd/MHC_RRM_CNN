from pandas import Series
from utils import load_data
from iedb import mhc_datasets

def get_matching_data(peplearn,
                      second,
                      allele='HLA',
                      peptide_col='peptide',
                      mhc_col='mhc'):
    '''
    Given a dataset with peptides and alleles from the peplearn database, and
    a second dataset containing peptides and alleles, returns a Series with
    the indexes of all rows in the seconds dataset which contain matching
    alleles and peptides to the peplearn dataset.

    Arguments -------
    peplearn: dataset from peplearn database
    second: other dataset of alleles and peptides
    peptide_col: optional, name of column in second with peptides defaults to
    'peptide'
    mhc_col: optional, name of column in second with alleles defaults to
    'mhc'
    '''
    indices = []
    for index, allele in enumerate(peplearn[2]):
        if allele.split('-')[0] == 'HLA':
            indices.append(index)
    seqs = []
    for index in indices:
        seq = (second[peptide_col] == peplearn[0][index])
        allele = (second[mhc_col] == peplearn[2][index])
        seqs += list(second[seq & allele].index)
    return Series(seqs)

def generate_tts_datasets(data_loc, arr_cols=None, names=None):
    '''
    Given the location of a dataset, splits the dataset into train and
    test files based on whether each peptide is train or test
    in the peplearn dataset. Saves the splits of the dataset to csvs.

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