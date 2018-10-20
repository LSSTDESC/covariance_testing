# covariance_testing

# Instructions: 
Please generate a 3x2pt Non-Gaussian covariance matrix for the parameters below. Except for area this corresponds to a DES Y1 covariance matrix assuming a circular survey geometry scaled to 5000 deg^2. Note that in the covariance released as part of the DES Y1 release the actual DES footprint is assumed in order to compute the correct noise terms. Here we neglect the specific DES footprint.  

# Goal: 
Code comparison at the precision required for an approximate DES Y3 analysis

## Cosmology
Omega_m : 0.2837.  
Omega_v : 0.7163

sigma_8 : 0.795431

n_spec : 0.96859

w0 : -1.0

wa : 0.0

omb : 0.062

h0 : 0.8433

## Angular Binning tmin,tmax in arcminutes
source_tomobins : 4

lens_tomobins : 5

shear_REDSHIFT_FILE : data/cosmolike_cov_Y3ish/source.nz

clustering_REDSHIFT_FILE : data/cosmolike_cov_Y3ish//lens.nz

tmin : 2.5

tmax : 250.0

ntheta : 20

## Survey Parameters

area : 5000.0

sigma_e : 0.39441348451

source_n_gal : 1.496,1.5189,1.5949,0.7949

lens_n_gal : 0.0134,0.0343,0.0505,0.0301,0.0089

lens_tomogbias : 1.44,1.70,1.698,1.997,2.058

