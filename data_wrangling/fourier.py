import numpy as np
import matplotlib.pyplot as plt

EIIP_DICT = {'L':0.0000,'I':0.0000,'N':0.0036,'N':0.0036,'G':0.0050,'V':0.0057,
             'E':0.0058,'P':0.0198,'H':0.0242,'K':0.0371,'A':0.0373,'Y':0.0516,
             'W':0.0548,'Q':0.0761,'M':0.0823,'S':0.0829,'C':0.0829,'T':0.0941,
             'F':0.0954,'R':0.0956,'D':0.1263}
             
HEMOGLOBIN_ALPHA = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYLLLLLL"

# Takes in an array of strings (the amino acid sequences) and returns a
# numerical array that corresponds to its EIIP values.
def eiip(seq):
    # If sequence length is uneven, add one from start
    if (len(seq) % 2) != 0.0:
        seq = seq + seq[0:1]
    x = []
    N = len(seq)
    # Creates data set of amino acid eiip values
    for s in range( N ):
        x.append(EIIP_DICT[seq[s]])
    return x


# Takes in an array of numerical EIIP values and returns the numerical
# information to create a Fourier transform of the paramaterized frequencies.
def fourierOf(eiipSeq):
    x = eiipSeq
    N = len(eiipSeq)
    n = np.reshape(np.linspace(1, int(N/2), int(N/2)), [1,int(N/2)])
    m = np.reshape(np.linspace(1,int(N), int(N)), [1, int(N)])
    v = np.matmul(np.transpose(m), n)
    # cos(x) portion
    cosx = np.cos(2 * np.pi/N * v)
    # sin(x) portion
    sinx = np.sin(2 * np.pi/N * v)
    # Gen eq: y(x) = R(Xre*cos(x) + Xim*sin(x))
    # Real part of X
    Xre = np.reshape(np.matmul(x, cosx), [1, int(N/2)])
    # Imaginary part of X
    Xim = np.reshape(np.matmul(x, sinx), [1, int(N/2)])
    # Amplitude: sqrt(Xre^2 + Xim^2)
    R = np.reshape(np.sqrt(Xre**2 + Xim**2),[1,int(N/2)])
    # Phase Shift: 
    phi = np.reshape(np.arctan(Xim/Xre), [1, int(N/2)])
    # Scaling
    scaled = 100*(R-min(min(R)))/(max(max(R))-min(min(R)))
    return n, N, Xre, Xim, R, phi, m, x, scaled