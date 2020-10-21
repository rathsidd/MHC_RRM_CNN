'''
File contains methods for working with MHC data.
'''

import numpy as np
import torch
import pandas as pd
import matplotlib.pyplot as plt


def convert_array(array):
    '''
    Takes in a list cast to a string and returns a numpy or pytorch array.

    Arguments ----
    array: a string containing an array or list of numbers
    '''
    array = array.strip('[]').split(', ')
    array = map(float, array)
    return np.array(list(array), dtype=np.float)


def load_data(csv, auto_convert=False, arr_cols=None, names=None):
    '''
    Takes path to csv containing MHC data, returns the values as a pd.DataFrame.

    Arguements ----
    csv: path to csv containing MHC data as a string.
    auto_convert: (optional) boolean, if True method will automatically convert
    strings it thinks are actually arrays to arrays.
    arr_cols: (optional) a list of names of columns containing arrays casted
    to strings which will be converted back to arrays.
    names: (optional) list of names to save converted columns as.
    '''
    data = pd.read_csv(csv)
    if arr_cols != None:
        for i, col in enumerate(arr_cols):
            if col not in list(data.columns):
                print('Error, column ' + arr_cols[i] + ' is not in dataset, continuing...')
            elif type(data[col][0]) != str:
                print('Error, column ' + arr_cols[i] + ' is not a string, continuing...')
            elif data[col][0][0] != '[':
                print('Error, column ' + arr_cols[i] + ' is not an array, continuing...')
            else:
                name = col
                if names != None:
                    name = names[i]
                data[name] = data[col].apply(convert_array)
    if auto_convert:
        for i, col in enumerate(data.columns):
            if (arr_cols != None and col not in arr_cols) or arr_cols == None:
                if type(data[col][0]) == str and data[col][0][0] == '[':
                    name = '' + col + '_as_array'
                    data[name] = data[col].apply(convert_array)
    return data


def plot_folds(ys, preds):
    '''
    Takes in lists of equal amounts of y_test values and corresponding model
    predictions and plots them as historgrams.

    Arguments ----
    ys: list of y_test arrays
    preds: list of prediction arrays
    '''
    fig, ax = plt.subplots(len(ys), figsize=(20,30))
    bins = np.linspace(-20, 10, 100)
    for i, axis in enumerate(ax):
	    axis.hist(preds[i],bins, alpha=0.5, label='preds')
	    axis.hist(ys[i],bins, alpha=0.5, label='ys')
	    axis.legend(loc='upper right')
	    axis.title.set_text(('Predictions vs Actual pIC50 values for fold ' + str(i + 1)))
    plt.savefig('linreg_results.png')