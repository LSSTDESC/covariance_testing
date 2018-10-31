import sys
infile = sys.argv[1]
outfile = sys.argv[2]

import twopoint
import numpy as np
data = twopoint.TwoPointFile.from_fits(infile)

covmat = data.covmat
print covmat
L = np.linalg.cholesky(covmat)
r = np.random.randn(len(L))

noise = np.dot(L,r)

i = 0
for s in data.spectra:
    n = len(s)
    s.value += noise[i:i+n]
    i+=n
data.to_fits(outfile)