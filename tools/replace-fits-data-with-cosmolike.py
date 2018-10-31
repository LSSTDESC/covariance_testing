#replace-cosmosis-datavector.py
import twopoint
import numpy as np
import sys


xip_order = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]
xim_order = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]
gammat_order = [(1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4),(5,1), (5,2),(5,3),(5,4),]
wtheta_order = [(1,1),(2,2),(3,3),(4,4),(5,5)]

cosmolike_order=[('xip', xip_order), ('xim', xim_order), ('gammat', gammat_order), ('wtheta', wtheta_order)]


def replace(fits_filename, cosmolike_filename, output_filename):
	fits_data = twopoint.TwoPointFile.from_fits(fits_filename)
	cosmolike_data = np.loadtxt(cosmolike_filename)
	if cosmolike_data.ndim==2:
		cosmolike_data = cosmolike_data[:,1]

	cosmosis_length = sum(len(s) for s in fits_data.spectra)
	cosmolike_length = len(cosmolike_data)
	if cosmolike_length!=cosmosis_length:
		raise ValueError("Data vectors are different lengths: cosmosis={} and cosmolike={}".format(cosmosis_length,cosmolike_length))


	total_bin_pairs = sum(len(x[1]) for x in cosmolike_order)
	angle_bins = cosmosis_length // total_bin_pairs

	i = 0
	#reordering the cosmosis data vector to match the cosmolike one
	for kind, order_list in cosmolike_order:
		order = []
		spectrum = fits_data.get_spectrum(kind)
		for d1,d2 in order_list:
			#Note the swap in order here to match the cosmolike one
			if kind.startswith("xi"):
				d1,d2=d2,d1
			w = np.where((spectrum.bin1==d1) & (spectrum.bin2==d2))[0]
			spectrum.value[w] = cosmolike_data[i:i+angle_bins]
			i+=angle_bins

	fits_data.to_fits(output_filename)



if __name__ == '__main__':
	fits_filename = sys.argv[1]
	cosmolike_filename = sys.argv[2]
	output_filename = sys.argv[3]
	replace(fits_filename, cosmolike_filename, output_filename)
