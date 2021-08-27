#####FIT OF VMB
## RUN:
## python3 fit_VMB.py 


import scipy
from scipy.optimize import curve_fit
import numpy as np

#Loading Average VMB
V = np.loadtxt('AVE_VMB')

# Initializing guess
p0 = np.zeros(2)
p0[0] = 7
p0[1] = 0.001

## VMB FIT
param, p_e = scipy.optimize.curve_fit(lambda t,p1,p2: p1*np.log(1+p2*t), V[:,0], V[:,1] , p0)

print("a", param[0],"b", param[1])

