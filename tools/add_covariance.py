"""
This is a template for an example script showing how you 
could add a covariance matrix to an existing file.
"""
import sys

twopoint_filename = sys.argv[1]
covmat_filename = sys.argv[2]
new_filename = sys.argv[3]


import twopoint
import numpy as np



print "Loading 2pt data in {}".format(twopoint_filename)
data = twopoint.TwoPointFile.from_fits(twopoint_filename, covmat_name=None)



print "Loading text covariance file from {}".format(covmat_filename)
covmat_data = np.loadtxt(covmat_filename).T
bin1 = covmat_data[0].astype(int)
bin2 = covmat_data[1].astype(int)
if len(covmat_data)==3:
	covmat_vector = covmat_data[2]
	print "Three columns - using third as cov"
elif len(covmat_data)==10:
	covmat_vector = covmat_data[8] + covmat_data[9]
	print "Ten columns - using sum of last two as cov"
else:
	print "Unknown format!"
	sys.exit(1)

print "Now you need to make sure that the ordering of the covariance matrix"
print "matches the ordering in the text file!"
print "Expecting this order:"
print "#Name Bin1 Bin2 Angle"
n = 0
for s in data.spectra:
	for b1, b2, ang in zip(s.bin1, s.bin2, s.angle):
		# print n, s.name, b1, b2, ang
		n+=1
print "Total data vector length: ", n

covmat = np.zeros((n,n))
for (b1,b2,c) in zip(bin1,bin2,covmat_vector):
	covmat[b1,b2] = c
	covmat[b2,b1] = c

cosmolike_xi_order = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]

#reordering the cosmosis data vector to match the cosmolike one
for kind in ['xip', 'xim']:
	order = []
	spectrum = data.get_spectrum(kind)
	for d1,d2 in cosmolike_xi_order:
		#Note the swap in order here to match the cosmolike one
		w = np.where((spectrum.bin1==d2) & (spectrum.bin2==d1))[0]
		order.append(w)
	order = np.concatenate(order)
	spectrum.apply_mask(order)

names = [s.name for s in data.spectra]
lengths = [len(s) for s in data.spectra]
n = sum(lengths)

assert covmat.shape==(n,n)

data.covmat_info = twopoint.CovarianceMatrixInfo("COVMAT", names, lengths, covmat)
data.to_fits(new_filename, clobber=True)
