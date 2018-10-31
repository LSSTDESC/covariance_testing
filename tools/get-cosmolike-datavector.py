"""
This is a template for an example script showing how you 
could add a covariance matrix to an existing file.
"""
import twopoint
import numpy as np
import sys


def setup(options):
	from cosmosis.datablock import option_section
	filename = options[option_section, "filename"]
	outfile = options[option_section, "outfile"]
	return filename, outfile

def execute(block, config):
	filename, outfile = config
	convert(filename, outfile)
	return 0


def convert(twopoint_filename, output_filename):
	data = twopoint.TwoPointFile.from_fits(twopoint_filename, covmat_name=None)


	xip_order = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]
	xim_order = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]
	gammat_order = [(1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4),(5,1),(5,2),(5,3),(5,4),]
	wtheta_order = [(1,1),(2,2),(3,3),(4,4),(5,5)]


	cosmolike_order=[('xip', xip_order), ('xim', xim_order), ('gammat', gammat_order), ('wtheta', wtheta_order)]
	        


	#reordering the cosmosis data vector to match the cosmolike one
	for kind, order_list in cosmolike_order:
		order = []
		spectrum = data.get_spectrum(kind)
		for d1,d2 in order_list:
			#Note the swap in order here to match the cosmolike one
			if kind.startswith("xi"):
				d1,d2=d2,d1
			w = np.where((spectrum.bin1==d1) & (spectrum.bin2==d2))[0]
			order.append(w)
		order = np.concatenate(order)
		spectrum.apply_mask(order)
	f = open(output_filename, "w")
	i=0
	for s in data.spectra:
		for v in s.value:
			f.write("%d  %e\n"%(i,v))
			i+=1
	f.close()




if __name__ == '__main__':
	twopoint_filename = sys.argv[1]
	output_filename = sys.argv[2]
	convert(twopoint_filename, output_filename)
