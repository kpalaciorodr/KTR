#####FIT OF VMB
## RUN:
## python3 fit_VMB.py 


import scipy
from scipy.optimize import curve_fit
import numpy as np

#Loading Average VMB
T=np.loadtxt('VMB')[:,0]             # TIME FOR ALL TRAJECTORIES (CONCATENATED)
V=np.loadtxt('VMB')[:,1]             # VMB FOR ALL TRAJECTORIES (CONCATENATED)

unique_T = np.unique(T)
unique_V = np.empty_like(unique_T)
for n_t, t in enumerate(unique_T):
    unique_V[n_t] = V[np.where(T == t)].mean()

# Initializing guess
p0 = np.zeros(2)
p0[0] = 7
p0[1] = 0.001

## VMB FIT
param, p_e = scipy.optimize.curve_fit(lambda t,p1,p2: p1*np.log(1+p2*t), unique_T, unique_V, p0)

print("a", param[0],"b", param[1])

