"""
np.conv():
    
pep = [peptide eiip array]
prot = [protein eiip array]
np.conv(pep, prot)

Dataframe Format:
    
index, peptide (eiip), protein_allele (eiip), conv (length), meas
"""

import os
import sys
import pandas as pd
import numpy as np
from numpy import asarray
from numpy import save
from numba import jit, cuda
from timeit import default_timer as timer

rawdata = pd.read_csv("non-convolved-CNN_data.csv")

def RRM(seq):
    if (len(seq)%2) != 0.0:
        seq = seq + seq[:1]
    N = len(seq)
    lib =  {"L":0.0000,
		"I":0.0000,
		"N":0.0036,
		"G":0.0050,
		"V":0.0057,
		"E":0.0058,
		"P":0.0198,
		"H":0.0242,
		"K":0.0371,
		"A":0.0373,
		"Y":0.0516,
		"W":0.0548,
		"Q":0.0761,
		"M":0.0823,
		"S":0.0829,
		"C":0.0829,
		"T":0.0941,
		"F":0.0946,
		"R":0.0959,
		"D":0.1263}
    x = []
    for s in range(N):
        x.append(lib[seq[s]])
    return x

count = 0
for seq in rawdata["peptide"]:
    rawdata.loc[count, "peptide_EIIP"] = str(RRM(rawdata.loc[count, "peptide"]))
    rawdata.loc[count, "protein_EIIP"] = str(RRM(rawdata.loc[count, "protein"]))
    count= count + 1

rawdata.to_csv('CNN_with_EIIP.csv', index = False)

sys.path.append('.')


#x = mhc_slot_dataset(allele=allele, directory=directory)

#df = {'Allele': train_alleles,
#      'Peptide Seq': train_seq,
 #     'IC50 Value': train_pIC50s,
 #     'Protein Seq': mhc_proteins}

#df = pd.DataFrame(df)
#df.to_csv('CNN_Data.csv', index=False)


