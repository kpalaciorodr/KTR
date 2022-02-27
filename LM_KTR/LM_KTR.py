#####NUMERICAL FIT OF VMB AND LIKELIHOOD OPTIMIZATION
## RUN:
## python3 LM_KTR.py <t_total> <VMB_file> <FINALTIMES_file> <sanity_checks(True, False), optional>

### INPUTS AND PARAMETERS
#t_total         = Total simulation time
#VMB_file        = File containing average maximum bias over all the simulations
#FINALTIMES_file = File containing vector with final simulation time (for non-events) or jump time (for events)
#sanity_checks   = Run sanity checks, writing additional files

import numpy as np
from sys import argv
import glob
from scipy import interpolate, optimize, integrate
import warnings
import multiprocessing as mp
from functools import partial
import pandas as pd

warnings.filterwarnings('ignore')

# Definition of functions

# Value of the loglikelihood as a function of gamma
# find mean_t directly from solving the dL/dmean_t=0 equation
def calculate_log_l(gamma, event, t, spline):

    # Run cum_hazard in parallel, change the 4 for something else if you want more cores
    p =  mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()
    
    log_hazard = calculate_log_hazard(gamma, t, spline)
    
    mean_t = cum_hazard.sum() / event.sum()
    log_l = -event.sum() * np.log(mean_t) + log_hazard[event].sum() - (1 / mean_t) * cum_hazard.sum()

    return -log_l

def calculate_cum_hazard(gamma, spline, t):

    int_Veff = integrate.quad(lambda x: np.exp(gamma * spline(x)), 0, t)[0]
    return int_Veff

def calculate_log_hazard(gamma, t, spline):

    Veff = spline(t)
    return gamma * Veff

# Compute CDF for our model using fitted parameters
def cumulative_density_fct(gamma, mean_t, spline, t):
    # Run cum_hazard in parallel, change the 4 for something else if you want more cores
    p = mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()
    return 1 - np.exp(-cum_hazard / mean_t)


# Perform one sample  Kolmogorov-Smirnov test
def ks_1samp_test(gamma, mean_t, spline, event, t):
    # Only take times with events
    N = np.sum(event)
    x = np.sort(t)[:N]
    cdfvals = cumulative_density_fct(gamma, mean_t, spline, x)
    Dplus = (np.arange(1.0, N + 1) / N - cdfvals).max()
    Dminus = (cdfvals - np.arange(0.0, N) / N).max()
    D = np.max([Dplus, Dminus])
    prob = stats.distributions.kstwo.sf(D, N)
    prob = np.clip(prob, 0, 1)
    return D, prob  # {"statistic": D, "pvalue": prob}

# Estimate the likelihood - objetive function
def ll_objective_function(spline,t,event):
    gamma_range = np.linspace(0, 1, 100)
    obj = np.zeros_like(gamma_range)
    for n, b in enumerate(gamma_range):
        obj[n] = calculate_log_l(b, event, t, spline)
    return gamma_range, obj

def main():
    # Load data and parameters
    #t=np.loadtxt('FINALTIMES')      
    t_total=int(argv[1])            # TOTAL SIMULATION TIME
    VMB=str(argv[2])                # NAME AVERAGE VMB FILE
    FINALTIMES=str(argv[3])         # NAME FINAL TIMES FILE
    t=np.loadtxt(FINALTIMES)[:,0]
    unique_T = np.loadtxt(VMB)[:,0] # TIME TO FIT VMB(t)
    unique_V = np.loadtxt(VMB)[:,1] # AVERAGE VMB TO FIT VMB(t)

    ## SPLINE FIT VMB(t)
    spline = interpolate.UnivariateSpline(unique_T, unique_V, s=0, ext=3)
    event = t < t_total

    ## LIKELIHOOD OPTIMIZATION
    opt = optimize.minimize_scalar(calculate_log_l, bounds=(0.0, 1.0), method='bounded', args=(event, t, spline))
    gamma = opt.x

    p =  mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()

    mean_t = cum_hazard.sum() / event.sum()
    result = np.array([1/mean_t, gamma])

    ## EXTRACTED GAMMA AND K
    print("gamma", result[1],"k", result[0])

    if not run_sanity_checks:
        return
    
    # Save VMB fit
    t_range=np.linspace(0,np.max(t),1000)
    V_range=spline(t_range)
    np.savetxt('spline_{}'.format(VMB),np.column_stack((t_range,V_range)))

    # Save empirical CDF
    counts = np.sort(t)
    norm = t.size
    ecdf = np.arange(1, norm + 1) / norm
    np.savetxt('ecdf_{}'.format(VMB), np.column_stack((counts, ecdf)))

    # Save fitted CDF from likelihood
    t_range = np.linspace(np.min(t), np.max(t), 150)
    cdfvals = cumulative_density_fct(gamma, mean_t, spline, t_range)
    np.savetxt('tcdf_{}'.format(VMB), np.column_stack((t_range, cdfvals)))

    # Save likelihood - objetive function
    np.savetxt('obj_likelihood_{}'.format(VMB),np.column_stack(ll_objective_function(spline,t,event)))

if __name__ == "__main__":

    run_sanity_checks = True if len(argv) >= 5 and argv[4].lower() == "true" else False
    main()

