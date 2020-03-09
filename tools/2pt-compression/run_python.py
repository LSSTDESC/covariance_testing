import os
import string
import numpy as np
from astropy.io import fits
import sys, tarfile
import string
from scipy import interpolate
import subprocess
import ConfigParser
import cosmosis.datablock

class Weight(object):
	def __init__(self, datafile, delta, **kwargs):
		self.param = ['omega_m', 'h0', 'log1e10As', 'omega_b', 'n_s', 'omnuh2', \
					'b1', 'b2', 'b3', 'b4', 'b5', \
					'm1', 'm2', 'm3', 'm4', \
					'A', 'alpha',\
					'bias_1', 'bias_2', 'bias_3', 'bias_4', \
					'ias_1', 'ias_2', 'ias_3', 'ias_4', 'ias_5']
		# The values coincide with what you are setting in the values.ini file
		self.value = [3.317044e-01, 7.273940e-01, 3.247298, 4.948368e-02, 9.782277e-01, 5.049024e-03, \
					1.45, 1.55, 1.65, 1.8, 2.0, \
					1.255649e-02, 1.353598e-02, 1.998341e-03, 1.713709e-02, \
					9.112135e-01, 1.889537, \
					3.213998e-03, -2.047676e-02, 7.529327e-03, -1.664299e-02, \
					0.0, 0.0, 0.0, 0.0, 0.0]
		# The fiducial parameter values
		self.fiducial_pv = {self.param[i]:self.value[i] for i in range(0,len(self.param))}
		self.cosmological_parameters = ['omega_m', 'h0', 'log1e10As', 'omega_b', 'n_s', 'omnuh2', \
					'w', 'massive_nu', 'massless_nu', 'omega_k', 'tau', 'wa']
		self.bin_bias = ['b1','b2','b3','b4','b5']
		self.shear_calibration_parameters = ['m1','m2','m3','m4']
		self.intrinsic_alignment_parameters = ['A','alpha']
		self.wl_photoz_errors = ['bias_1','bias_2','bias_3','bias_4']
		self.lens_photoz_errors = ['ias_1','ias_2','ias_3','ias_4','ias_5']
		self.script_dir = os.path.dirname(__file__)
		self.output_dir = 'derivative-values/data_vector/'
		self.datafile = datafile
		self.delta = delta
		self.params = kwargs # a dictionary
		for key in kwargs:
		   setattr(self, key, kwargs[key])
		self.theta, self.data, self.cov, self.cov_inverse = self.get_data()
		self.weight =  self.get_weight()
		self.y_d = self.get_y()
		self.compressed_fidicual = self.get_y_t()

	def set_params(self,**kwargs):  #set param values in values-derivative.ini file
		values_dir = os.path.join(self.script_dir, "values-derivative.ini")
		config = ConfigParser.RawConfigParser()
		config.read(values_dir)
		for kw in kwargs:
			if kwargs[kw]==0.:
				delta = 0.0001
			else:
				delta = np.abs(self.delta*kwargs[kw])
			if kw in self.cosmological_parameters:
				config.set('cosmological_parameters','%s'%kw,'%e'\
					%(kwargs[kw]))
			elif kw in self.bin_bias:
				config.set('bin_bias','%s'%kw,'%f'\
					%(kwargs[kw]))
			elif kw in self.shear_calibration_parameters:
				config.set('shear_calibration_parameters','%s'%kw,'%f'\
					%(kwargs[kw]))
			elif kw in self.intrinsic_alignment_parameters:
				config.set('intrinsic_alignment_parameters','%s'%kw,'%f'\
					%(kwargs[kw]))
			elif kw in self.wl_photoz_errors:
				config.set('wl_photoz_errors','%s'%kw,'%f'\
					%(kwargs[kw]))
			elif kw in self.lens_photoz_errors:
				config.set('lens_photoz_errors','%s'%('b'+kw),'%f'\
					%(kwargs[kw]))
			with open(values_dir,'wb') as configfile:
				config.write(configfile)

	def run_pipeline(self):
		"""
		Input:: dictionary of cosmological parameters
		Output:: A list of (len_theory ,2) arrays, each array is two cosmosis
				theory vectors with perturbations p/m delta from fiducial values
		"""
		variations = {}   ## {key: parameter_name   value: perturbed theory}
		first = 0
		for kw in self.params:
			if self.fiducial_pv[kw]==0.:
				delta = 0.0001
			else:
				delta = np.abs(self.delta*self.params[kw])
			if first == 0:
				values = [self.params[kw]-delta,self.params[kw]+delta,self.params[kw]]
			else:
				values = [self.params[kw]-delta,self.params[kw]+delta]
			fiducial = {kw:self.params[kw]}
			theory_p = []
			for p in values:
				dictionary = {kw:p}
				self.set_params(**dictionary)
				os.system('cosmosis params_derivative.ini')
				theory_p.append(np.loadtxt(self.output_dir+'2pt_theory.txt'))
			theory_p = np.asarray(theory_p)
			print('theory shape:', np.shape(theory_p))
			if first == 0:
				self.fiducial_theory = theory_p.T[:,2]
			variations[kw] = theory_p
			first += 1
			self.reset(**fiducial)
		return variations

	def derivative(self):
		variations = self.run_pipeline()
		derivatives = {}
		for key in variations:
			if self.fiducial_pv[key]==0.:
				delta = 0.0001
			else:
				delta = np.abs(self.delta*self.params[key])
			deriv = (variations[key].T[:,1] - variations[key].T[:,0]) / (delta * 2.)
			derivatives[key] = deriv
		print('Finished derivative calculation')
		self.derivatives = derivatives
		return derivatives

	def get_weight(self):
		derivatives = self.derivative()
		weight = {}
		for key in derivatives:
			new_dev = derivatives[key]
			weight[key] = new_dev.dot(self.cov_inverse)
			weight[key]= (weight[key])/np.linalg.norm(weight[key])
		sorted_weight = sorted([[k,v] for k,v in weight.items()], key=lambda x: x[0])
		weight_array = np.asarray([x[1] for x in sorted_weight])
		return weight_array

	def get_y(self):
		y=np.dot(self.weight,self.data)
		return y

	def get_y_t(self):
		y_t = np.dot(self.weight,self.fiducial_theory)
		return y_t

	def get_data(self):
		"""
		Output: theta, data, covariance matrix and inverse covariance from the data file
		"""
		try:
			hdulist=fits.open(self.datafile)
		except ValueError:
			print("Could not open: ", self.datafile)
		theta = hdulist['xip'].data['ang'][:20]
		data = np.loadtxt(self.output_dir+'2pt_data.txt')
		cov = np.loadtxt(self.output_dir+'2pt_covariance.txt')
		cov_inv = np.loadtxt(self.output_dir+'2pt_inverse_covariance.txt')
		print('Got data from .fits file')
		for kw in self.params:
			print("We are calculating the derivatives for the following parameters:")
			print(kw)
		return theta, data, cov, cov_inv

	def reset(self,**kwargs):
		values_dir=os.path.join(self.script_dir, "values-derivative.ini")
		config= ConfigParser.RawConfigParser()
		config.read(values_dir)
		for kw in kwargs:
			if kw in self.cosmological_parameters:
				config.set('cosmological_parameters','%s'%kw,'%e'%kwargs[kw])
			elif kw in self.bin_bias:
				config.set('bin_bias','%s'%kw,'%f'%kwargs[kw])
			elif kw in self.shear_calibration_parameters:
				config.set('shear_calibration_parameters','%s'%kw,'%f'%kwargs[kw])
			elif kw in self.intrinsic_alignment_parameters:
				config.set('intrinsic_alignment_parameters','%s'%kw,'%f'%kwargs[kw])
			elif kw in self.wl_photoz_errors:
				config.set('wl_photoz_errors','%s'%kw,'%f'%kwargs[kw])
			elif kw in self.lens_photoz_errors:
				config.set('lens_photoz_errors','%s'%('b'+kw),'%f'%kwargs[kw])
			print("returned",'%s = '%kw,'%f'%kwargs[kw])
			with open(values_dir,'wb') as configfile:
				config.write(configfile)
		return 0

def values(filename, delta, **kwargs):
	trial1 = Weight(filename,delta,**kwargs)
	weight = trial1.weight
	y_d = trial1.y_d
	cov = trial1.cov
	theta = trial1.theta
	derivatives = trial1.derivatives
	fiducial_theory = trial1.fiducial_theory
	compressed_fidicual = trial1.compressed_fidicual
	cov_c = weight.dot(cov).dot(weight.T)
	return derivatives, cov_c, weight, y_d, theta, fiducial_theory, compressed_fidicual
