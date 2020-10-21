# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:59:27 2020

@author: Andrew
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ast

rawdata = pd.read_csv("CNN_with_EIIP.csv")


count = 0
for seq in rawdata["conv_length"]:
    x = ast.literal_eval(rawdata.loc[count, "peptide_EIIP"])
    y = ast.literal_eval(rawdata.loc[count, "protein_EIIP"])
    rawdata.loc[count, "conv_length"] = str(np.convolve(x, y, mode = 'same').tolist())
    count= count + 1

rawdata.to_csv('CNN_with_EIIP.csv', index = False)



# =============================================================================
# x = [0.0829, 0.0829, 0.1263, 0.0373, 0.0946, 0.0516, 0.0198, 0.0946, 0.0516, 0.0829]
# y = []
# 
# convolution_same = np.convolve(x, y, mode = 'same')
# axis = np.linspace(0, len(y), len(y))
# plt.plot(axis, convolution_same)
# plt.xlabel("domain")
# plt.ylabel("same conv.")
# plt.show()
# 
# convolution_full = np.convolve(x, y, mode = 'full')
# axis = np.linspace(0, len(convolution_full), len(convolution_full))
# plt.plot(axis, convolution_full)
# plt.xlabel("domain")
# plt.ylabel("full conv.")
# plt.show()
# =============================================================================
