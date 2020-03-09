from run_python import Weight, values
import numpy as np
import pickle
import os

""" FIDUCIAL VALUES FOR REFERENCE:
omega_m = 2.932322e-01
h0 = 7.337380e-01
omega_b = 4.839119e-02
n_s = 9.839522e-01
log1e10As = 3.247298
omnuh2 = 5.068874e-03
m1 = 1.205911e-02
m2 = 1.384025e-02
m3 = 2.800220e-03
m4 = 1.776998e-02
a = 8.360370e-01
alpha = 2.100269
bias_1 = 2.548395e-03
bias_2 = -2.008459e-02
bias_3 = 8.146909e-03
bias_4 = -1.609595e-02
"""

os.system('cosmosis params_derivative.ini')

save_dir = '../Parameters'
derivatives, cov_c, weight, y, theta, fiducial_theory, compressed_fiducial= \
values('2pt_NG.fits', \
								0.005, omega_m = 2.932322e-01, h0 = 7.337380e-01, omega_b = 4.839119e-02, n_s = 9.839522e-01, \
								omnuh2 = 5.068874e-03, log1e10As = 3.247298, \
								m1 = 1.205911e-02, m2 = 1.384025e-02, m3 = 2.800220e-03, m4 = 1.776998e-02, \
								A = 8.360370e-01, alpha = 2.100269, \
								bias_1 = 2.548395e-03, bias_2 = -2.008459e-02, bias_3 = 8.146909e-03, bias_4 = -1.609595e-02)

derivative = derivatives.values()
derivative = np.asarray(derivative)

print(weight)
print(y)
print(cov_c)

try:
    os.stat(save_dir)
except:
    os.mkdir(save_dir)

np.savetxt(save_dir+'/weight.txt',weight)
np.savetxt(save_dir+'/data_c.txt',y)
np.savetxt(save_dir+'/cov_c.txt',cov_c)
np.savetxt(save_dir+'/theta.txt',theta)
np.savetxt(save_dir+'/derivatives.txt',derivative)
np.savetxt(save_dir+'/fiducial_theory.txt',fiducial_theory)
np.savetxt(save_dir+'/compressed_fiducial.txt',compressed_fiducial)
