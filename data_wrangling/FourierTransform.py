# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:17:17 2020

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt

eiip_dict = {'L':0.0000,'I':0.0000,'N':0.0036,'N':0.0036,'G':0.0050,'V':0.0057,
             'E':0.0058,'P':0.0198,'H':0.0242,'K':0.0371,'A':0.0373,'Y':0.0516,
             'W':0.0548,'Q':0.0761,'M':0.0823,'S':0.0829,'C':0.0829,'T':0.0941,
             'F':0.0954,'R':0.0956,'D':0.1263}

# Human hemoglobin alpha
seq1 = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYLLLLLL"
# Human hemoglobin beta
seq2 = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
# Chimpanzee hemoglobin beta
seq3 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRLLLL"
# Rat hemoglobin beta
seq4 = "MVHLTDAEKAAVNGLWGKVNPDDVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVINAFNDGLKHLDNLKGTFAHLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKEFTPCAQAAFQKVVAGVASALAHKY"

def fourierOf(seq):
    # If sequence length is uneven, add one from start
    if (len(seq)%2) != 0.0:
        seq = seq + seq[:1]
    x = []
    N = len(seq)
    # Creates data set of amino acid eiip values
    for s in range(N):
        x.append(eiip_dict[seq[s]])
    n = np.reshape(np.linspace(1,int(N/2),int(N/2)),[1,int(N/2)])
    m = np.reshape(np.linspace(1,int(N),int(N)),[1,int(N)])
    v = np.matmul(np.transpose(m),n)
    # cos(x) portion
    cosx = np.cos(2*np.pi/N*v)
    # sin(x) portion
    sinx = np.sin(2*np.pi/N*v)
    #Gen eq: y(x) = R(Xre*cos(x) + Xim*sin(x))
    #Real part of X
    Xre = np.reshape(np.matmul(x,cosx), [1,int(N/2)])
    #Imaginary part of X
    Xim = np.reshape(np.matmul(x,sinx), [1,int(N/2)])
    #Amplitude: sqrt(Xre^2 + Xim^2)
    R = np.reshape(np.sqrt(Xre**2 + Xim**2),[1,int(N/2)])
    #Phase Shift: 
    phi = np.reshape(np.arctan(Xim/Xre),[1,int(N/2)])
    #...
    scaled = 100*(R-min(min(R)))/(max(max(R))-min(min(R)))
    return n, N, Xre, Xim, R, phi, m, x, scaled




result1 = fourierOf(seq1)

result2 = fourierOf(seq2)

#result3 = fourierOf(seq3)

#result4 = fourierOf(seq4)

#S = result1[2]*result2[2]*result3[2]*result4[2] + result1[3]*result2[3]*result3[3]*result4[3]
S = result1[2]*result2[2] + result1[3]*result2[3]
#S = result1[2] + result1[3]

S = 100*(S-min(min(S)))/(max(max(S))-min(min(S)))

# SEQUENCE
plt.scatter(result1[6],np.reshape(result1[7],[1,146]))
plt.show()
plt.scatter(result2[6],np.reshape(result2[7],[1,146]))
plt.show()

# FREQ
plt.plot(np.reshape(result1[0]/result1[1],[int(result1[1]/2),]),np.reshape(result1[8],[int(result1[1]/2),]))
plt.axis([0, 0.5, 0, 100])
plt.show()

plt.plot(np.reshape(result2[0]/result2[1],[int(result2[1]/2),]),np.reshape(result2[8],[int(result2[1]/2),]))
plt.axis([0, 0.5, 0, 100])
plt.show()

# CROSS SPECTRUM/SPECTRA
plt.plot(np.reshape(result1[0]/result1[1],[int(result1[1]/2),]),np.reshape(100*S/max(max(S)),[int(result1[1]/2),]))
plt.axis([0, 0.5, 0, 100])
plt.show()
