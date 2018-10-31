#generate-posterior-sample.py

import sys
infile = sys.argv[1]
outfile = sys.argv[2]

import numpy as np
data = np.loadtxt(infile)
nvaried = None
nsample = None
for line in open(infile):
    if not line.startswith('#'):
        continue
    if line.startswith('#n_varied'):
        nvaried = int(line.split('=')[-1])
    if line.startswith('#nsample'):
        nsample = int(line.split('=')[-1])

if nvaried is None:
    raise ValueError("Input chain does not have a line #nvaried=... near the start that multinest chains should have")
if nsample is None:
    raise ValueError("Input chain does not have a line #nsample=... near the end that multinest chains should have")


data = data[-nsample:]
weight = data[:,-1]
weight = weight/weight.max()
p = np.random.uniform(size=len(weight))
accept = weight>p
data = data[accept,:nvaried]


np.savetxt(outfile, data)
