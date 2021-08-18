#####KTR method to extract k and gamma from biased simulations
## RUN:
## python3 cal_KTR.py <n> <a> <b>

### INPUTS AND PARAMETERS
#t = vector with simulation jump times for each event
#a = input from log fit of VMB = a*log(1+b*t)
#b = input from log fit of VMB = a*log(1+b*t)
#p0 = initial guesses for optimizing the fit 
#n = total number of simulations


import scipy
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
from sys import argv
from numpy.random import random

#Loading Jumptimes
t=np.loadtxt('JUMPTIMES')

mintime=min(t)
maxtime=max(t)

#Reading n, a1 and b1 params from log fit of VMB
n=int(argv[1])
a1=argv[2]
b1=argv[3]

a=float(a1)
b=float(b1)

# Initializing guess
p0=np.zeros(2)

gamma0=0.6 ## initial gamma that enables extracting a converged fit

p0[0]=gamma0*a+1.0
p0[1]=0.0000000001

## Calculating empirical CDF
logbins = np.linspace(mintime,maxtime,n)

counts, bin_edges = np.histogram(t, bins=logbins)

ecdf = np.cumsum(counts)

##For printing ECDF
#for i in range(0,len(ecdf)):
#       print('CDF',logbins[i],ecdf[i]/float(n))

## ECDF FIT
param,p_e= scipy.optimize.curve_fit(lambda t,p1,p2: 1-np.exp(-(pow(1+b*t,p1)-1.0)*p2), logbins[0:-1],  ecdf/float(n), p0, bounds=(0, [p0[0]+3,np.inf]))

### CREATING THEORETICAL CDF
cdf = lambda x:  1-np.exp(-(pow(1.0+b*x,param[0])-1.0)*param[1])

binscdf = np.linspace(0,0.95*maxtime,2000)
cdfbin= np.zeros(2000)

for i in range(0,2000):
      cdfbin[i]=cdf(binscdf[i])

theoretical_times=[]

for i in range(0,10000):
       c = random( )
       n = 0
       for j in range(0,2000):
             if cdfbin[j]>c and n==0:
                   n=1
                   theoretical_times.append(binscdf[j])

## CALCULATING KS TEST
s,p_val = scipy.stats.ks_2samp(t,theoretical_times)

## EXTRACTED GAMMA, K AND P-VALUE
print("gamma", (param[0]-1.0)/a,"k", param[1]*param[0]*b, "pval",p_val)
