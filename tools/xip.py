# The goal is to read in Niall's new format and input into
# the pipeline so the relevant n(z)'s and 2-point functions
# are there
import numpy as np
import matplotlib.pyplot as plt
import fitsio
plt.style.use('SVA1StyleSheet.mplstyle')

###
# ang index  a: 0:14
# zindex z: 11:0 12: 1 13: 2 14: 3 22: 4 23: 5 24: 6 33: 7 34: 8 44: 9 

nz=10
ndesa = 20
na = 15 # from BJ
thmin = 8.

########################
## Read in the xip block from the fits file
########################
filename = '/Users/sdodelso/des-mpp/data_vectors/y1/2pt_NG_mcal_1110.fits'
tp_fits = fitsio.FITS(filename)
xip = tp_fits['xip']
theta = xip['ang'][:ndesa]
thc = np.delete(theta,range(ndesa-na))

ccl = tp_fits[1]
ccl_xip=ccl[:ndesa*nz,:ndesa*nz]
print ccl_xip[199,199]
cclcut=np.zeros((nz*na,nz*na))


for z1 in range(nz):
	for a1 in range(na):
		i1 = a1+z1*na
		## The index in the original matrix is skewed because that matrix contains
		## angles that have been cut. So 
		i1ccl = (a1+ndesa-na) + z1*ndesa
		for z2 in range(nz):
			for a2 in range(na):
				i2 = a2 + z2*na
				## The index in the original matrix is skewed because that matrix contains
				## angles that have been cut. So 
				i2ccl = (a2+ndesa-na) + z2*ndesa
				cclcut[i1,i2]=ccl_xip[i1ccl,i2ccl]

########################
## Read in the xip block from the BJ file
########################
bjcov = np.loadtxt('thps_cov_des_matrix.dat').T
## Throw away xim which is every second block of na
bjxip=np.zeros((nz*na,nz*na))
for z1 in range(nz):
	for a1 in range(na):
		i1 = a1+z1*na
		## The index in the original matrix is skewed because that matrix contains
		## xim 
		i1bjcov = (a1) + z1*na*2
		for z2 in range(nz):
			for a2 in range(na):
				i2 = a2 + z2*na
				i2bjcov = (a2) + z2*na*2
				bjxip[i1,i2] = bjcov[i1bjcov,i2bjcov]

print np.max(bjxip),np.min(bjxip)
plt.scatter(cclcut,bjxip)
plt.plot([1.e-14,1.e-11],[1.e-14,1.e-11])
plt.axis([1.e-14,1.e-10,1.e-14,1.e-10])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$C_{CL}$')
plt.ylabel('$C_{BJ}$')
plt.savefig('xipscatter.png',bbox_inches='tight')
				
